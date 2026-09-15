require 'bibtex'

module Jekyll
  # Generates one standalone page per publication at /publications/<key>/,
  # carrying Google Scholar's Highwire Press citation_* meta tags (see
  # _includes/metadata.html). Scholar's spec expects these on a page about
  # exactly one article; the main /publications/ page lists all of them, so
  # it can't correctly carry citation tags for any single paper.
  #
  # Each page is a legitimate standalone bibliographic record (title,
  # authors, venue, DOI, abstract, PDF link) rather than a bare stub, and
  # bib_alternative.html links each entry's title to its page so these are
  # reachable through normal site navigation, not just the sitemap.
  class PublicationPageGenerator < Generator
    priority :normal

    def generate(site)
      scholar = site.config['scholar'] || {}
      source_dir = scholar['source'] || '_bibliography'
      bib_file = scholar['bibliography'] || 'papers.bib'
      bib_path = File.join(site.source, source_dir.sub(%r{\A/}, ''), bib_file)
      return unless File.exist?(bib_path)

      BibTeX.open(bib_path).each do |entry|
        next unless entry.is_a?(BibTeX::Entry)

        key = entry.key.to_s
        title = clean(entry[:title])
        next if title.empty?

        authors = entry[:author] ? entry[:author].to_a : []
        citation_authors = authors.map { |n| [n.last, n.first].reject { |s| s.to_s.empty? }.join(', ') }
        authors_display = authors.map { |n| [n.first, n.last].reject { |s| s.to_s.empty? }.join(' ') }.join(', ')

        is_conference = entry.type == :inproceedings
        # Some journal/booktitle values in the .bib already end in a period
        # (e.g. "arXiv."); strip it here so the layout's own "." doesn't
        # double up, and citation_journal_title stays clean for Scholar.
        venue = clean(is_conference ? entry[:booktitle] : entry[:journal]).sub(/\.+\z/, '')

        year = entry[:year] ? entry[:year].to_s : nil
        month = month_number(entry[:month])
        citation_date = month ? "#{year}/#{month}" : year

        doi = entry[:doi] ? entry[:doi].to_s.strip : nil
        doi = nil if doi && doi.empty?

        pdf = entry[:pdf] ? entry[:pdf].to_s.strip : nil
        pdf = nil if pdf && pdf.empty?
        citation_pdf_url = direct_pdf_url?(pdf) ? pdf : nil

        page = PageWithoutAFile.new(site, site.source, "publications/#{key}", 'index.html')
        page.content = ''
        page.data.merge!(
          'layout' => 'citation',
          'title' => title,
          'permalink' => "/publications/#{key}/",
          'description' => "#{title} — #{authors_display}",
          'pub_key' => key,
          'pub_type' => is_conference ? 'Conference paper' : 'Journal article',
          'authors_display' => authors_display,
          'venue' => venue,
          'year' => year,
          'volume' => entry[:volume] ? entry[:volume].to_s : nil,
          'pages' => entry[:pages] ? entry[:pages].to_s : nil,
          'doi' => doi,
          'pdf_url' => pdf,
          'abstract_text' => clean(entry[:abstract]),
          'citation_title' => title,
          'citation_authors' => citation_authors,
          'citation_type' => is_conference ? 'conference' : 'journal',
          'citation_venue' => venue,
          'citation_date' => citation_date,
          'citation_doi' => doi,
          'citation_pdf_url' => citation_pdf_url
        )
        site.pages << page
      end
    end

    private

    MONTH_NUMBERS = {
      'jan' => '01', 'feb' => '02', 'mar' => '03', 'apr' => '04',
      'may' => '05', 'jun' => '06', 'jul' => '07', 'aug' => '08',
      'sep' => '09', 'oct' => '10', 'nov' => '11', 'dec' => '12'
    }.freeze

    # bibtex-ruby normalizes a bare numeric `month = 9` into the BibTeX
    # month macro (:sep), so entry[:month].to_s returns "sep", not "9" -
    # map either form to a zero-padded number for citation_publication_date.
    def month_number(value)
      return nil if value.nil?

      raw = value.to_s.strip.downcase
      return nil if raw.empty?
      return raw.rjust(2, '0') if raw =~ /\A\d{1,2}\z/

      MONTH_NUMBERS[raw[0, 3]]
    end

    # Strips the curly braces BibTeX uses for case-preservation, and
    # collapses whitespace (abstracts in this .bib file are hard-wrapped).
    def clean(value)
      return '' if value.nil?
      value.to_s.gsub(/[{}]/, '').gsub(/\s+/, ' ').strip
    end

    def direct_pdf_url?(url)
      return false if url.nil? || url.empty?
      return false if url.include?('github.com') && url.include?('/blob/')
      return true if url =~ /\.pdf(\?.*)?\z/i
      return true if url =~ %r{/pdf/?\z}i
      return true if url =~ %r{arxiv\.org/pdf/}i

      false
    end
  end
end
