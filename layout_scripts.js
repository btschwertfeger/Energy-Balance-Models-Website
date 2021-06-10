/*
################################################################################################################################
########################################################################################################
##################################################################################
################### SIMPLE EBM
*/
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

/*
################################################################################################################################
########################################################################################################
##################################################################################
################### COMPLEX EBM
*/

let complex_ebm_submitbtn = document.getElementById("complex_ebm_submitbtn");
complex_ebm_submitbtn.onclick = function () {
    window.doTheComplexThing({
        D: document.getElementById('complex_input_D').value,
        S1: document.getElementById('complex_input_S1').value,
        A: document.getElementById('complex_input_A').value,
        B: document.getElementById('complex_input_B').value,
        cw: document.getElementById('complex_input_cw').value,
        S0: document.getElementById('complex_input_S0').value,
        S2: document.getElementById('complex_input_S2').value,
        a0: document.getElementById('complex_input_a0').value,
        a2: document.getElementById('complex_input_a2').value,
        ai: document.getElementById('complex_input_ai').value,
        Fb: document.getElementById('complex_input_Fb').value,
        k: document.getElementById('complex_input_k').value,
        Lf: document.getElementById('complex_input_Lf').value,
        cg: document.getElementById('complex_input_cg').value,
        tau: document.getElementById('complex_input_tau').value,
        winter: document.getElementById('complex_input_winter').value,
        summer: document.getElementById('complex_input_summer').value,
        years: document.getElementById('complex_input_years').value
    });
};

let complex_variables = ["D", "S1", "A", "B", "cw", "S0", "S2", "a0", "a2", "ai", "Fb", "k", "Lf", "cg", "tau", "winter", "summer", "years"];
let complex_ebm_inputFields = document.getElementsByName('complex-ebm_input');

let complex_ebm_resettbtn = document.getElementById("complex_ebm_resetbtn");
complex_ebm_resettbtn.onclick = function () {
    complex_variables.forEach((element, index) => {
        let default_value = window.default_complex_ebm_input[element];
        document.getElementById("complex_input_" + element).value = default_value;
    });
}
// window.doTheComplexThing(window.default_complex_ebm_input);