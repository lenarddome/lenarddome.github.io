---
layout: page
title: software
permalink: /software/
description: I build software for research and data analysis. All of them are open-source and free to use. You will never be asked to pay for any of these. Scientist should not be charged for using tools to push what we know further <i class="far fa-smile"></i>
nav: true
nav_order: 3
---

<script data-name="BMC-Widget" data-cfasync="false" src="https://cdnjs.buymeacoffee.com/1.0.0/widget.prod.min.js" data-id="lenarddome" data-description="Support me on Buy me a coffee!" data-message="" data-color="#FFDD00" data-position="Right" data-x_margin="18" data-y_margin="18"></script>

<!-- Place this tag where you want the button to render. -->
<a class="github-button" href="https://github.com/lenarddome" data-size="large" data-show-count="true" aria-label="Follow @lenarddome on GitHub">Follow @lenarddome</a>

{% if site.data.repositories.github_users %}
<div class="repositories d-flex flex-wrap flex-md-row flex-column justify-content-between align-items-center">
  {% for user in site.data.repositories.github_users %}
    {% include repository/repo_user.html username=user %}
  {% endfor %}
</div>
{% endif %}

---

## GitHub Software Repositories

{% if site.data.repositories.github_repos %}
<div class="repositories d-flex flex-wrap flex-md-row flex-column align-items-left">
  {% for repo in site.data.repositories.github_repos %}
    {% include repository/repo.html repository=repo %}
  {% endfor %}
</div>
{% endif %}

