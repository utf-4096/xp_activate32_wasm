import XPA from './xpa.js';

const webAssemblyIsSupported = ('WebAssembly' in window);
const xpa = new XPA();
let keyIn, keyOut, btnGenerate, elmErrorContainer, elmError;

function showError(error) {
    if (!error) {
        elmErrorContainer.style.display = 'none';
        return;
    }

    elmErrorContainer.style.display = '';
    elmError.textContent = error.toString();
}

function generateKey() {
    try {
        keyOut.value = xpa.generate(keyIn.value);
        showError(null);
    } catch(err) {
        showError(err.toString());
    }
}

document.addEventListener('DOMContentLoaded', () => {
    /* get elements */
    keyIn = document.querySelector('#key-in');
    keyOut = document.querySelector('#key-out');
    btnGenerate = document.querySelector('#generate');
    elmErrorContainer = document.querySelector('#error-container');
    elmError = document.querySelector('#error');

    if (!webAssemblyIsSupported) {
        showError('WebAssembly is required to use this program');
        return;
    }

    /* events */
    btnGenerate.addEventListener('click', generateKey);
    keyIn.addEventListener('keyup', ({ key }) => {
        if (key === 'Enter') {
            generateKey();
        }
    });
});

if (webAssemblyIsSupported) {
    await xpa.boot();
}
