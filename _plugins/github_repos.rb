require 'net/http'
require 'json'
require 'fileutils'

module Jekyll
  # Fetches live repo metadata (stars, description, language, license) from
  # the GitHub API for each repo in _data/repositories.yml, so the software
  # page can render real cards instead of depending on a third-party
  # card-image service (github-readme-stats.vercel.app, which was returning
  # 503s in production).
  #
  # Exposes site.data['github_repos'] - the same list, in the same order,
  # each entry augmented with 'stars', 'description', 'language', 'license',
  # 'html_url', 'homepage', and (when a matching _software/<name>.markdown
  # exists) 'detail_url' pointing at its page on this site.
  #
  # Network calls happen once per build. A local cache
  # (.jekyll-cache/github_repos.json) is read as a fallback whenever the API
  # call fails - rate limited, offline, GitHub down - so a flaky network
  # never breaks the build; it just serves the last known-good data. If
  # there's no cache either (first build, offline), falls back to the bare
  # owner/repo slug so the page still renders something.
  class GithubReposGenerator < Generator
    priority :normal
    CACHE_PATH = File.expand_path('.jekyll-cache/github_repos.json', __dir__ + '/../')
    TIMEOUT = 6

    def generate(site)
      repos = site.data.dig('repositories', 'github_repos') || []
      return if repos.empty?

      cache = load_cache
      software_slugs = software_detail_urls(site)

      fetched = repos.map do |slug|
        data = fetch(slug) || cache[slug] || fallback(slug)
        data['detail_url'] = software_slugs[data['name'].to_s.downcase]
        data
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

    def fetch(slug)
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
        'slug' => slug,
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
      Jekyll.logger.warn 'GithubRepos:', "#{slug} fetch failed (#{e.class}: #{e.message}), falling back"
      nil
    end

    def fallback(slug)
      owner, name = slug.split('/')
      {
        'slug' => slug,
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
