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

Needs the Ruby version in `.ruby-version` (e.g. `chruby 3.4.10`).

```bash
bundle install
bundle exec jekyll serve
```

## Deploy

Pushing to `source` runs `.github/workflows/deploy.yml`, which builds the site
and publishes it to GitHub Pages. Dependabot opens monthly PRs for gem and
Actions updates.

Fonts, icons and the code-highlighting theme are self-hosted under
`assets/vendor/` and `assets/css/highlight/`; MathJax only loads on pages with
`math: true` in their front matter.
