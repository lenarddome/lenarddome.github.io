---
layout: page
permalink: /software/
title: software
description: Some selected software and code.
nav: true
nav_order: 1
---

I am keen on developing scientific software. In Psychology and probably Cognitive Science, this type of endeavour is often underappreciated and underfunded. However, it is a crucial part of the scientific process. I am a strong advocate for open-source science and I am involved in a variety of open-source science projects. Here are some of the software projects I am involved in:

---

## GitHub Repositories

{% if site.data.github_repos %}
<div class="repo-grid">
  {% for repo in site.data.github_repos %}
    {% include repository/repo_card.html repository=repo %}
  {% endfor %}
</div>
{% endif %}
