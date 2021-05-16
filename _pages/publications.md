---
layout: page
permalink: /publications/
title: publications
description: 
years: [2021, 2019]
nav: true
---

For a list that includes other resources and supplementary
materials, visit the github repository: [publications](https://github.com/lenarddome/publications)

<div class="publications">

{% for y in page.years %}
  <h3 class="year">{{y}}</h3>
  {% bibliography -f papers -q @*[year={{y}}]* %}
{% endfor %}

</div>
