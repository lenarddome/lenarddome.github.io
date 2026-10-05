# lenarddome.com

Source for [lenarddome.com](https://lenarddome.com), a Jekyll site built on the
[al-folio](https://github.com/alshedivat/al-folio) theme.

## Layout

- `_pages/` — top-level pages (about, publications, software, news, videos, resume)
- `_bibliography/papers.bib` — publications, rendered by jekyll-scholar
- `_data/` — coauthors, keyword colours, and the repositories shown on `/software/`
- `_news/`, `_posts/` — news items and blog posts
- `_plugins/` — GitHub repo cards, per-publication citation pages, keyword histogram

## Local development

```bash
bundle install
bundle exec jekyll serve
```

## Deploy

Work happens on the `source` branch. `bin/deploy` builds the site and pushes it
to `gh-pages`, which GitHub Pages serves.
