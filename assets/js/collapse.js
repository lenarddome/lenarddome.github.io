// Expands and collapses `.collapse` elements, replacing Bootstrap's collapse
// plugin: a `[data-toggle="collapse"][data-target="#id"]` button toggles #id,
// animating its height via the `.collapsing` class in _sass/_framework.scss.
(function () {
    function transitionMs(element) {
        var duration = parseFloat(getComputedStyle(element).transitionDuration) || 0;
        return duration * 1000;
    }

    // Calls `done` once the height transition ends (or straight away when
    // transitions are off, e.g. for visitors who prefer reduced motion).
    function afterTransition(element, done) {
        var ms = transitionMs(element);
        if (!ms) {
            done();
            return;
        }
        var finished = false;
        function finish() {
            if (finished) return;
            finished = true;
            element.removeEventListener('transitionend', finish);
            done();
        }
        element.addEventListener('transitionend', finish);
        setTimeout(finish, ms + 50);
    }

    function show(element) {
        element.classList.remove('collapse');
        element.classList.add('collapsing');
        element.style.height = '0px';
        element.offsetHeight; // force a reflow so the height change animates
        element.style.height = element.scrollHeight + 'px';
        afterTransition(element, function () {
            element.classList.remove('collapsing');
            element.classList.add('collapse', 'show');
            element.style.height = '';
        });
    }

    function hide(element) {
        element.style.height = element.getBoundingClientRect().height + 'px';
        element.offsetHeight;
        element.classList.add('collapsing');
        element.classList.remove('collapse', 'show');
        element.style.height = '';
        afterTransition(element, function () {
            element.classList.remove('collapsing');
            element.classList.add('collapse');
        });
    }

    document.addEventListener('click', function (event) {
        var trigger = event.target.closest('[data-toggle="collapse"]');
        if (!trigger) return;
        var target = document.querySelector(trigger.getAttribute('data-target'));
        if (!target || target.classList.contains('collapsing')) return;
        event.preventDefault();

        var opening = !target.classList.contains('show');
        if (opening) {
            show(target);
        } else {
            hide(target);
        }
        document.querySelectorAll('[data-toggle="collapse"][data-target="' + trigger.getAttribute('data-target') + '"]')
            .forEach(function (button) {
                button.classList.toggle('collapsed', !opening);
                button.setAttribute('aria-expanded', opening ? 'true' : 'false');
            });
    });
})();
