---
layout: page
permalink: /publications/
title: publications
description:
years: [2023, 2021, 2019]
nav: true
nav_order: 1
---
<!-- _pages/publications.md -->
<div class="publications">

<hr>
<br>

{%- for y in page.years %}
  {% bibliography -f {{ site.scholar.bibliography }} -q @*[year={{y}}]* %}
{% endfor %}

</div>
