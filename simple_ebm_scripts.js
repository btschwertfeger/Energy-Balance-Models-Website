/* GENERAL VARIABLES*/
let tebm_sliders = document.getElementsByName("tebm_slider");
let tebm_variabels = ['D', 'A', 'B', 'cw', 'S0', 'S2', 'a0', 'a2', 'ai', 'F', 'gamma'];


/* FINAL YEAR CHART */
/* ADD BUTTON  */
let tebm_addbtn = document.getElementById('add_tebm_graph');
tebm_addbtn.onclick = function () {
    window.added_tebm_graphs += 1;
}
/* RESET BUTTON  */
let tebm_resetbtn = document.getElementById('tebm_resetbtn');
tebm_resetbtn.onclick = function () {
    window.plot_default_tebm_chart();
    tebm_sliders.forEach((element, index) => {
        let default_value = window.default_TEBM_input[tebm_variabels[index]];
        document.getElementById(element.id).value = default_value;
        document.getElementById(tebm_variabels[index] + '_sliderAmount').innerHTML = default_value;
    });
    window.added_tebm_graphs = 0;
    document.getElementById('xlatitudes-container').value = 'all';
}
/* SLIDER  */
for (let entry = 0; entry < tebm_sliders.length; entry++) {
    tebm_sliders[entry].onchange = function () {
        window.doTEBM({
            D: document.getElementById('D_slide').value,
            A: document.getElementById('A_slide').value,
            B: document.getElementById('B_slide').value,
            cw: document.getElementById('cw_slide').value,
            S0: document.getElementById('S0_slide').value,
            S2: document.getElementById('S2_slide').value,
            a0: document.getElementById('a0_slide').value,
            a2: document.getElementById('a2_slide').value,
            ai: document.getElementById('ai_slide').value,
            F: document.getElementById('F_slide').value,
            gamma: document.getElementById('gamma_slide').value
        })
        document.getElementById(tebm_sliders[entry].id + 'rAmount').innerHTML = document.getElementById(tebm_sliders[entry].id).value;
    }
}