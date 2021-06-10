import {
    array,
    first
} from '/Users/benjamin/js_modules/HTMLCollection/HTMLCollection.js';

document.getElementsByClassName("child").array().forEach((item, index) => {
    item.innerHTML = "TEST" + index;
});