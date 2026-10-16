// Copy buttons: <button data-copy-target="id"> copies that element's text.
document.addEventListener('click', function (event) {
    const button = event.target.closest('[data-copy-target]');
    if (!button) return;
    const target = document.getElementById(button.dataset.copyTarget);
    navigator.clipboard.writeText(target.textContent).then(function () {
        const icon = button.querySelector('i');
        button.classList.add('copied');
        icon.className = 'fa-solid fa-check';
        button.title = 'Copied';
        setTimeout(function () {
            button.classList.remove('copied');
            icon.className = 'fa-regular fa-copy';
            button.title = 'Copy BibTeX';
        }, 1500);
    }, function () {
        // Clipboard access refused: select the text so it can be copied by hand.
        window.getSelection().selectAllChildren(target);
    });
});
