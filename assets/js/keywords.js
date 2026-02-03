(function () {
    document.addEventListener('DOMContentLoaded', function () {
        var pills = Array.prototype.slice.call(document.querySelectorAll('.keyword-pill'));
        if (!pills.length) {
            return;
        }

        var lockedKey = null;

        function applyHighlight(key) {
            pills.forEach(function (pill) {
                var matches = pill.dataset.keyword === key;
                pill.classList.toggle('is-highlighted', matches);
                pill.classList.toggle('is-muted', !matches);
            });
        }

        function clearHighlight() {
            lockedKey = null;
            pills.forEach(function (pill) {
                pill.classList.remove('is-highlighted');
                pill.classList.remove('is-muted');
            });
        }

        function handleEnter(event) {
            if (lockedKey) {
                return;
            }
            var key = event.currentTarget.dataset.keyword;
            if (!key) {
                return;
            }
            applyHighlight(key);
        }

        function handleLeave() {
            if (lockedKey) {
                return;
            }
            clearHighlight();
        }

        function handleFocus(event) {
            if (lockedKey) {
                return;
            }
            var key = event.currentTarget.dataset.keyword;
            if (!key) {
                return;
            }
            applyHighlight(key);
        }

        function handleBlur() {
            if (lockedKey) {
                return;
            }
            clearHighlight();
        }

        function handleClick(event) {
            var key = event.currentTarget.dataset.keyword;
            if (!key) {
                return;
            }
            if (lockedKey === key) {
                clearHighlight();
            } else {
                lockedKey = key;
                applyHighlight(key);
            }
        }

        pills.forEach(function (pill) {
            pill.addEventListener('mouseenter', handleEnter);
            pill.addEventListener('mouseleave', handleLeave);
            pill.addEventListener('focus', handleFocus);
            pill.addEventListener('blur', handleBlur);
            pill.addEventListener('click', handleClick);
        });

        document.addEventListener('click', function (event) {
            if (!event.target.closest('.keyword-pill')) {
                if (lockedKey) {
                    clearHighlight();
                }
            }
        });

        document.addEventListener('keydown', function (event) {
            if (event.key === 'Escape' && lockedKey) {
                clearHighlight();
            }
        });
    });
})();
