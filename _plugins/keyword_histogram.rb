require 'bibtex'

module Jekyll
  # Builds a keyword -> frequency histogram from the same .bib file
  # jekyll-scholar renders, so it can never drift out of sync with the
  # publication list. Exposes `site.data.keyword_histogram`, an array of
  # { keyword, slug, count, entries: [{ key, title, year }] }, sorted by
  # count desc then alphabetically. `slug` matches the `keyword_slug` used
  # for `data-keyword` in _layouts/bib_alternative.html (both go through
  # Jekyll::Utils.slugify(..., mode: 'raw')), so the histogram can drive
  # the same highlighting/filtering as the keyword pills.
  class KeywordHistogramGenerator < Generator
    priority :normal

    def generate(site)
      scholar = site.config['scholar'] || {}
      source_dir = scholar['source'] || '_bibliography'
      bib_file = scholar['bibliography'] || 'papers.bib'
      bib_path = File.join(site.source, source_dir.sub(%r{\A/}, ''), bib_file)

      unless File.exist?(bib_path)
        site.data['keyword_histogram'] = []
        return
      end

      counts = Hash.new(0)
      entries_by_keyword = Hash.new { |h, k| h[k] = [] }

      BibTeX.open(bib_path).each do |entry|
        next unless entry.is_a?(BibTeX::Entry)
        next unless entry[:keywords]

        title = entry[:title] ? entry[:title].to_s.gsub(/[{}]/, '') : entry.key.to_s
        year = entry[:year] ? entry[:year].to_s : nil

        entry[:keywords].to_s.split(',').each do |raw|
          clean = raw.strip
          next if clean.empty?

          key = clean.downcase
          counts[key] += 1
          entries_by_keyword[key] << { 'key' => entry.key.to_s, 'title' => title, 'year' => year }
        end
      end

      site.data['keyword_histogram'] = counts.map do |key, count|
        {
          'keyword' => key,
          'slug' => Utils.slugify(key, mode: 'raw'),
          'count' => count,
          'entries' => entries_by_keyword[key]
        }
      end.sort_by { |h| [-h['count'], h['keyword']] }
    end
  end
end
