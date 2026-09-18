require 'net/http'
require 'json'
require 'fileutils'

module Jekyll
  # Fetches live repo metadata (stars, description, language, license) from
  # the GitHub API, and - where _data/repositories.yml declares a published
  # package (`cran:`/`pypi:`) - a download count from that registry, for
  # each repo in _data/repositories.yml. Replaces a third-party card-image
  # service (github-readme-stats.vercel.app) that was returning 503s in
  # production.
  #
  # Package registry names are declared explicitly per repo rather than
  # guessed from the GitHub repo name: a PyPI project named "cpm" exists and
  # is a completely unrelated package, so name-guessing would have silently
  # attributed someone else's downloads to this project.
  #
  # Exposes site.data['github_repos'] - one entry per repo, in the same
  # order, with 'stars', 'description', 'language', 'license', 'html_url',
  # 'homepage', 'downloads', 'downloads_source' ("CRAN"/"PyPI", or nil when
  # no registry is declared or the count wasn't available), and (when a
  # matching _software/<name>.markdown exists) 'detail_url'.
  #
  # Network calls happen once per build. A local cache
  # (.jekyll-cache/github_repos.json) is read as a fallback whenever a call
  # fails - rate limited, offline, service down - so a flaky network never
  # breaks the build; it just serves the last known-good value for whichever
  # part failed. GitHub metadata falls back to a bare owner/repo link if
  # there's no cache either; a download count with no cache is simply
  # omitted rather than shown as a misleading zero.
  class GithubReposGenerator < Generator
    priority :normal
    CACHE_PATH = File.expand_path('.jekyll-cache/github_repos.json', __dir__ + '/../')
    TIMEOUT = 6

    def generate(site)
      entries = site.data.dig('repositories', 'github_repos') || []
      return if entries.empty?

      cache = load_cache
      software_slugs = software_detail_urls(site)

      fetched = entries.map do |entry|
        slug = entry['repo'] || entry.to_s
        cached = cache[slug] || {}

        github = fetch_github(slug) || extract_github_fields(cached) || fallback_github(slug)
        downloads = fetch_downloads(entry)
        downloads['downloads'] ||= cached['downloads']
        downloads['downloads_formatted'] ||= cached['downloads_formatted']
        downloads['downloads_source'] ||= cached['downloads_source']

        github.merge(downloads).merge(
          'slug' => slug,
          'detail_url' => software_slugs[github['name'].to_s.downcase]
        )
      end

      save_cache(fetched)
      site.data['github_repos'] = fetched
    end

    private

    # Maps a repo name (e.g. "psp") to its /software/psp/ page, when a
    # matching entry exists in the software collection.
    def software_detail_urls(site)
      docs = site.collections['software'] ? site.collections['software'].docs : []
      docs.each_with_object({}) do |doc, map|
        map[doc.basename_without_ext.downcase] = doc.url
      end
    end

    def fetch_github(slug)
      uri = URI("https://api.github.com/repos/#{slug}")
      http = Net::HTTP.new(uri.host, uri.port)
      http.use_ssl = true
      http.open_timeout = TIMEOUT
      http.read_timeout = TIMEOUT

      request = Net::HTTP::Get.new(uri)
      request['Accept'] = 'application/vnd.github+json'
      request['User-Agent'] = 'lenarddome.github.io-build'
      request['Authorization'] = "Bearer #{ENV['GITHUB_TOKEN']}" if ENV['GITHUB_TOKEN']

      response = http.request(request)
      return nil unless response.is_a?(Net::HTTPSuccess)

      json = JSON.parse(response.body)
      {
        'name' => json['name'],
        'owner' => slug.split('/').first,
        'description' => json['description'],
        'stars' => json['stargazers_count'],
        'language' => json['language'],
        'license' => json.dig('license', 'spdx_id'),
        'html_url' => json['html_url'],
        'homepage' => (json['homepage'].to_s.empty? ? nil : json['homepage'])
      }
    rescue StandardError => e
      Jekyll.logger.warn 'GithubRepos:', "#{slug} GitHub fetch failed (#{e.class}: #{e.message}), falling back"
      nil
    end

    GITHUB_FIELDS = %w[name owner description stars language license html_url homepage].freeze

    # The cache stores each repo's fields flattened together (see
    # save_cache); pull just the GitHub-sourced subset back out as a
    # fallback when a live fetch fails. Returns nil for an empty/missing
    # cache entry so the caller falls through to fallback_github instead of
    # a hash of all-nil fields.
    def extract_github_fields(cached)
      return nil if cached.nil? || cached.empty?

      cached.slice(*GITHUB_FIELDS)
    end

    def fallback_github(slug)
      owner, name = slug.split('/')
      {
        'name' => name,
        'owner' => owner,
        'description' => nil,
        'stars' => nil,
        'language' => nil,
        'license' => nil,
        'html_url' => "https://github.com/#{slug}",
        'homepage' => nil
      }
    end

    EMPTY_DOWNLOADS = { 'downloads' => nil, 'downloads_formatted' => nil, 'downloads_source' => nil }.freeze

    def fetch_downloads(entry)
      return EMPTY_DOWNLOADS.dup unless entry.is_a?(Hash)

      if entry['cran']
        fetch_cran_downloads(entry['cran'])
      elsif entry['pypi']
        fetch_pypi_downloads(entry['pypi'])
      else
        EMPTY_DOWNLOADS.dup
      end
    rescue StandardError => e
      Jekyll.logger.warn 'GithubRepos:', "download count failed (#{e.class}: #{e.message}), falling back"
      EMPTY_DOWNLOADS.dup
    end

    # cranlogs.r-pkg.org's grand-total range gives an exact all-time count -
    # matches what CRAN's own download badges show.
    def fetch_cran_downloads(package)
      total = get_json("https://cranlogs.r-pkg.org/downloads/total/2000-01-01:2099-12-31/#{package}")
                &.first&.dig('downloads')
      { 'downloads' => total, 'downloads_formatted' => format_number(total), 'downloads_source' => 'CRAN' }
    end

    # PyPI has no free, unauthenticated, exact all-time count: pypistats.org
    # is exact but only keeps a rolling ~180-day window (which undercounts a
    # package's real lifetime total), and pepy.tech's exact API needs a key.
    # pepy.tech's public badge *is* all-time and needs no key, but only
    # returns an abbreviated figure ("3k", "1.2m") - shown as-is rather than
    # converted to a fake-precise number, since the badge doesn't expose the
    # rounding it used.
    def fetch_pypi_downloads(package)
      svg = get_text("https://static.pepy.tech/badge/#{package}")
      values = svg&.scan(%r{<text[^>]*>([^<]*)</text>})&.flatten&.uniq
      label = values&.last

      formatted = label && (label =~ /[km]\z/i ? "~#{label.upcase}" : label)

      {
        'downloads' => parse_abbreviated_count(label),
        'downloads_formatted' => formatted,
        'downloads_source' => 'PyPI'
      }
    end

    # "3k" -> 3000, "1.2m" -> 1200000, "179" -> 179; nil/unparseable -> nil.
    def parse_abbreviated_count(label)
      return nil unless label

      match = label.strip.downcase.match(/\A([\d.]+)([km]?)\z/)
      return nil unless match

      multiplier = { 'k' => 1_000, 'm' => 1_000_000 }.fetch(match[2], 1)
      (match[1].to_f * multiplier).round
    end

    def get_json(url)
      body = get_text(url)
      body && JSON.parse(body)
    end

    def get_text(url)
      uri = URI(url)
      http = Net::HTTP.new(uri.host, uri.port)
      http.use_ssl = true
      http.open_timeout = TIMEOUT
      http.read_timeout = TIMEOUT

      request = Net::HTTP::Get.new(uri)
      request['User-Agent'] = 'lenarddome.github.io-build'
      response = http.request(request)
      return nil unless response.is_a?(Net::HTTPSuccess)

      response.body
    end

    # Comma-groups a number for display (18175 -> "18,175"); nil stays nil.
    def format_number(number)
      return nil if number.nil?

      number.to_s.reverse.gsub(/(\d{3})(?=\d)/, '\1,').reverse
    end

    def load_cache
      return {} unless File.exist?(CACHE_PATH)

      JSON.parse(File.read(CACHE_PATH)).each_with_object({}) do |entry, map|
        map[entry['slug']] = entry
      end
    rescue StandardError
      {}
    end

    def save_cache(entries)
      FileUtils.mkdir_p(File.dirname(CACHE_PATH))
      File.write(CACHE_PATH, JSON.pretty_generate(entries))
    rescue StandardError => e
      Jekyll.logger.warn 'GithubRepos:', "could not write cache (#{e.message})"
    end
  end
end
