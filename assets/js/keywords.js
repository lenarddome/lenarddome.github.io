(function () {
    document.addEventListener('DOMContentLoaded', function () {
        var controls = Array.prototype.slice.call(
            document.querySelectorAll('.keyword-pill, .keyword-histogram-bar')
        );
        if (!controls.length) {
            return;
        }

        var entries = Array.prototype.slice.call(document.querySelectorAll('.bibliography > li'));
        var lockedKey = null;

        function applyHighlight(key) {
            controls.forEach(function (control) {
                var matches = control.dataset.keyword === key;
                control.classList.toggle('is-highlighted', matches);
                control.classList.toggle('is-muted', !matches);
            });
        }

        function clearHighlight() {
            controls.forEach(function (control) {
                control.classList.remove('is-highlighted');
                control.classList.remove('is-muted');
            });
        }

        // Histogram bars sit above the whole list, so locking one also
        // filters the bibliography down to entries carrying that keyword.
        // Individual keyword pills only highlight (an entry can't hide itself).
        function applyFilter(key) {
            entries.forEach(function (entry) {
                var hasMatch = !!entry.querySelector('.keyword-pill[data-keyword="' + key + '"]');
                entry.classList.toggle('is-hidden', !hasMatch);
            });
        }

        function clearFilter() {
            entries.forEach(function (entry) {
                entry.classList.remove('is-hidden');
            });
        }

        function setPressed(key) {
            controls.forEach(function (control) {
                if (!control.classList.contains('keyword-histogram-bar')) {
                    return;
                }
                control.setAttribute('aria-pressed', control.dataset.keyword === key ? 'true' : 'false');
            });
        }

        function lock(key) {
            lockedKey = key;
            applyHighlight(key);
            applyFilter(key);
            setPressed(key);
        }

        function unlock() {
            lockedKey = null;
            clearHighlight();
            clearFilter();
            setPressed(null);
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
            handleEnter(event);
        }

        function handleBlur() {
            handleLeave();
        }

        function handleClick(event) {
            var key = event.currentTarget.dataset.keyword;
            if (!key) {
                return;
            }
            if (lockedKey === key) {
                unlock();
            } else {
                lock(key);
            }
        }

        controls.forEach(function (control) {
            control.addEventListener('mouseenter', handleEnter);
            control.addEventListener('mouseleave', handleLeave);
            control.addEventListener('focus', handleFocus);
            control.addEventListener('blur', handleBlur);
            control.addEventListener('click', handleClick);
        });

        document.addEventListener('click', function (event) {
            if (!event.target.closest('.keyword-pill, .keyword-histogram-bar')) {
                if (lockedKey) {
                    unlock();
                }
            }
        });

        document.addEventListener('keydown', function (event) {
            if (event.key === 'Escape' && lockedKey) {
                unlock();
            }
        });
    });
})();
