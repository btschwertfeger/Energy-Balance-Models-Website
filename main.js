// const {
//     fractionDependencies,
//     sec
// } = require("mathjs");
const {
    typeOf
} = require("mathjs");
const math = require("mathjs");
const linear = require("/Users/benjamin/js_modules/gauss-jordan.js");

/*###################################################################################################################################################*/
/*###################################################################################################################################################*/
/*###################################################################################################################################################*/
/*###################################################################################################################################################*/
/*###################################################################################################################################################*/
/*###################################################################################################################################################*/
/*####################################################### FUNCTIONS AND SETTINGS ####################################################################*/

window.COLORS = ['blue', 'black', 'orange', 'green', 'Chocolate', 'Crimson', 'DarkCyan', 'DarkGoldenRod', 'Fuchsia', 'Indigo', 'yellow'];

// returns array from start to end with steps 
function range(start, end, step) {
    var range = [];
    var typeofStart = typeof start;
    var typeofEnd = typeof end;

    if (step === 0) {
        throw TypeError("Step cannot be zero.");
    }

    if (typeofStart == "undefined" || typeofEnd == "undefined") {
        throw TypeError("Must pass start and end arguments.");
    } else if (typeofStart != typeofEnd) {
        throw TypeError("Start and end arguments must be of same type.");
    }
    typeof step == "undefined" && (step = 1);

    if (end < start) {
        step = -step;
    }
    if (typeofStart == "number") {
        while (step > 0 ? end >= start : end <= start) {
            range.push(start);
            start += step;
        }
    } else if (typeofStart == "string") {
        if (start.length != 1 || end.length != 1) {
            throw TypeError("Only strings with one character are supported.");
        }
        start = start.charCodeAt(0);
        end = end.charCodeAt(0);
        while (step > 0 ? end >= start : end <= start) {
            range.push(String.fromCharCode(start));
            start += step;
        }
    } else {
        throw TypeError("Only string and number types are supported");
    }
    return range;
}

// check if value is true for array.some(isTrue)
function isTrue(value) {
    return (value != NaN && value != 0 && value != undefined) ? true : false;
}

function isLessThanZero(value) {
    return (value < 0) ? true : false;
}

// returns array from to with equal distance
function linespace(startValue, stopValue, cardinality) {
    var arr = [];
    var step = (stopValue - startValue) / (cardinality - 1);
    for (var i = 0; i < cardinality; i++) {
        arr.push(startValue + (step * i));
    }
    return arr;
}

function dynamicColors() {
    let r = Math.floor(Math.random() * 255);
    let g = Math.floor(Math.random() * 255);
    let b = Math.floor(Math.random() * 255);
    return "rgba(" + r + "," + g + "," + b + ", 0.5)";
}

function poolColors(a) {
    let pool = [];
    for (let i = 0; i < a; i++) {
        pool.push(dynamicColors());
    }
    return pool;
}

/*###################################################################################################################################################*/
/*###################################################################################################################################################*/
/*###################################################################################################################################################*/
/*###################################################################################################################################################*/
/*###################################################################################################################################################*/
/*###################################################################################################################################################*/
/*############################################### TEMPERATURE EBM ###################################################################################*/


function calcTEBM(D = 0.6, A = 193, B = 2.1, cw = 9.8, S0 = 420, S2 = 240, a0 = 0.7, a2 = 0.1, ai = 0.4, F = 0.0, gamma = 1) {
    F = parseFloat(F); //sonst interpretiert er das komischerweise als String
    D = parseFloat(D);
    let n = 50 // grid resolution (number of points between equator and pole)
    let nt = .5
    let dur = 100
    let dt = 1 / nt
    // Spatial Grid ---------------------------------------------------------
    let dx = 1.0 / n // grid box width
    let x = range(dx / 2, 1 + dx / 2, dx);
    x = x.map(function (entry) {
        return Math.round(entry * 100) / 100
    }) //native grid
    let xb = range(dx, 1, dx);
    xb = xb.map(function (entry) {
        return Math.round(entry * 100) / 100
    });

    // Diffusion Operator (WE15, Appendix A) -----------------------------------
    let lam = (xb.map(function (entry) {
        return 1 - Math.pow(entry, 2)
    })).map(function (entry) {
        return entry * D / Math.pow(dx, 2)
    })

    let L1 = [];
    L1.push(0);
    lam.map(function (entry) {
        L1.push(Math.round(entry * -1 * 100) / 100)
    });
    let L2 = [];
    lam.map(function (entry) {
        L2.push(Math.round(entry * -1 * 100) / 100)
    });
    L2.push(0);
    let L3 = new Array(L1.length);
    for (var i = 0; i < L3.length; i++) {
        L3[i] = L1[i] * -1 - L2[i];
    };

    let L3_diag = [];
    for (let i = 0; i < L3.length; i++) {
        L3_diag.push(new Array(L3.length))
    }

    for (let row = 0; row < L3_diag.length; row++) {
        for (let column = 0; column < L3_diag.length; column++) {
            if (row === column) {
                L3_diag[row][column] = L3[row];
            } else {
                L3_diag[row][column] = 1;
                L3_diag[row][column] = 0;
            }
        }
    }
    let L2_diag = [];
    for (let i = 0; i < L2.length; i++) {
        L2_diag.push(new Array(L2.length))
    }
    for (let row = 0; row < L2_diag.length; row++) {
        for (let column = 0; column < L2_diag.length; column++) {
            if (row + 1 === column) {
                L2_diag[row][column] = L2[row];
            } else {
                L2_diag[row][column] = 1;
                L2_diag[row][column] = 0;
            }
        }
    }
    let L1_diag = [];
    for (let i = 0; i < L1.length; i++) {
        L1_diag.push(new Array(L1.length))
    }
    for (let row = 0; row < L1_diag.length; row++) {
        for (let column = 0; column < L1_diag.length; column++) {
            if (row === column + 1) {
                L1_diag[row][column] = L1[row];
            } else {
                L1_diag[row][column] = 1;
                L1_diag[row][column] = 0;
            }
        }
    }

    let diffop = [];
    for (let i = 0; i < L1.length; i++) {
        diffop.push(new Array(L1.length))
    }
    for (let row = 0; row < L1.length; row++) {
        for (let column = 0; column < L1.length; column++) {
            diffop[row][column] = -L3_diag[row][column] - L2_diag[row][column] - L1_diag[row][column];
        }
    }

    let S = x.map(function (entry) {
        return S0 - S2 * Math.pow(entry, 2)
    })
    let aw = x.map(function (entry) {
        return a0 - a2 * Math.pow(entry, 2)
    })

    let T = new Array(x.length).fill(10)
    let allT = [];
    for (let i = 0; i < dur * nt; i++) {
        allT.push(new Array(n))
    }
    for (let row = 0; row < dur * nt; row++) {
        for (let column = 0; column < n; column++) {
            allT[row][column] = 1;
            allT[row][column] = 0;
        }
    }
    let t = linespace(0, dur, Math.round(dur * nt));

    let I = []
    for (let i = 0; i < n; i++) {
        I.push(new Array(n));
    }
    for (let row = 0; row < n; row++) {
        for (let column = 0; column < n; column++) {
            I[row][column] = 1;
            if (row != column) {
                I[row][column] = 0;
            }
        }
    }
    /* I+dt/cw*(B*I-diffop) */

    let BI = []
    for (let i = 0; i < I.length; i++) {
        BI.push(new Array(I.length));
    }
    for (let row = 0; row < I.length; row++) {
        for (let column = 0; column < I.length; column++) {
            if (row == column) {
                BI[row][column] = B;
            } else {
                BI[row][column] = 1;
                BI[row][column] = 0;
            }
        }
    }

    // I+dt/cw*(BI-diffop)
    let BIdiffop = []
    for (let i = 0; i < I.length; i++) {
        BIdiffop.push(new Array(I.length));
    }

    //BIdiffop = I+dt/cw*(BI-diffop)
    for (let row = 0; row < BIdiffop.length; row++) {
        for (let column = 0; column < BIdiffop.length; column++) {
            BIdiffop[row][column] = BI[row][column] - diffop[row][column];
            BIdiffop[row][column] = BIdiffop[row][column] * (dt / cw);
            if (row == column) {
                BIdiffop[row][column] += I[row][column];
            }
        }
    }
    let invMat = math.inv(BIdiffop);

    for (let i = 0; i < parseInt(dur * nt); i++) {
        //python: a = aw*(T>0)+ai*(T<0) # WE15, eq.4 
        let awT = aw.map(function (entry, index) {
            if (T[index] > 0) {
                return entry;
            } else {
                return 0;
            }
        })

        let aiT = T.map(function (entry, index) {
            if (entry < 0) {
                return ai;
            } else {
                return 0;
            }
        });

        let a = awT.map(function (entry, index) {
            return entry + aiT[index];
        })

        //python: C = a*S-A+F 
        let C = a.map(function (entry, index) {
            return entry * S[index] - (gamma * A) + F;
        })

        //python:  T0 = T+dt/cw*C
        let T0 = C.map(function (entry) {
            return dt / cw * entry
        }).map(function (entry, index) {
            return entry + T[index]
        })

        // Governing equation [cf. WE15, eq. (2)]:
        // T(n+1) = T(n) + dt*(dT(n+1)/dt), with c_w*dT/dt=(C-B*T+diffop*T)
        // -> T(n+1) = T(n) + dt/cw*[C-B*T(n+1)+diff_op*T(n+1)]
        // -> T(n+1) = inv[1+dt/cw*(1+B-diff_op)]*(T(n)+dt/cw*C)
        T = math.multiply(invMat, T0)
        allT[i] = allT[i].map(function (entry, index) {
            return T[index]
        })

    }

    return {
        'allT': allT,
        't': t,
        'T': T,
        'x': x
    }
}

/*#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##-#-#-#-#-#*/

/* SCALE TO HAVE LATITUDE DATA FROM -82 TP +82 */
function scaleTEBMResults(data) {
    let invx = data['x'].map(function (entry) {
        return entry
    })
    data['x'].map(function (entry, index) {
        invx.unshift(-entry)
    })

    let invT = data["T"]
    let tmp = data['T'].map(function (entry, index) {
        return data['T'][data['T'].length - index - 1];
    })
    tmp.map(function (entry, index) {
        invT.unshift(tmp[tmp.length - index - 1]);
    });

    data['x'] = invx;
    data['T'] = invT;
    return data;
}

window.last_TEBM_res = [];

window.doTEBM = function doTEBM(TEBM_input_obj) {
    let D = TEBM_input_obj['D'];
    let A = TEBM_input_obj['A'];
    let B = TEBM_input_obj['B'];
    let cw = TEBM_input_obj['cw'];
    let S0 = TEBM_input_obj['S0'];
    let S2 = TEBM_input_obj['S2'];
    let a0 = TEBM_input_obj['a0'];
    let a2 = TEBM_input_obj['a2'];
    let ai = TEBM_input_obj['ai'];
    let F = TEBM_input_obj['F'];
    let gamma = TEBM_input_obj['gamma'];

    let TEBM_res = calcTEBM(D, A, B, cw, S0, S2, a0, a2, ai, F, gamma);
    TEBM_res = scaleTEBMResults(TEBM_res);
    window.last_TEBM_res = TEBM_res;
    updateTEBM_charts(); //TEBM_res['T'], TEBM_res['allT'])
}

/*#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##-#-#-#-#-#*/

/* UPDATES CHART AFTER SLIDER EVENT */
function updateTEBM_charts() { //T, allT) {
    /* SIMPLE CHART */
    let T = window.last_TEBM_res['T'];
    let allT = window.last_TEBM_res['allT'];

    let chart = window.tebm_chart;
    let graphCount = chart.data.datasets.length;

    if (graphCount > (1 + window.added_tebm_graphs)) {
        chart.data.datasets.pop();
        graphCount--;
    }

    const LABEL = 'Temperature by input #' + (graphCount);
    chart.data.datasets.push({
        label: LABEL,
        data: T,
        fill: false,
        borderColor: window.COLORS[chart.data.datasets.length],
        pointRadius: 0
    });

    if (chart.data.datasets.length > 7) {
        chart.options.plugins.legend.position = 'left';
    }

    chart.update();

    /* MULTICHART */
    let tebm_chart_all = window.tebm_chart_all;

    // let randColors = poolColors(allT.length)
    let randColors = poolColors(window.activeLatitudes.length)

    let allTbyLat = new Array(window.activeLatitudes.length);
    for (let i = 0; i < allTbyLat.length; i++) {
        allTbyLat[i] = new Array(allT.length); // add time dimension
    }

    for (let year = 0, idx = 0; year < allT.length; year++) {
        for (let lat = 0; lat < allT[year].length; lat++) {
            if (window.activeLatitudes.includes(lat)) {
                allTbyLat[idx][year] = allT[year][lat];
                idx++;
                idx %= window.activeLatitudes.length;
            }
        }
    }

    // let allTbyLat = new Array(allT[0].length)
    // for (let i = 0; i < allTbyLat.length; i++) {
    //     allTbyLat[i] = new Array(window.activeLatitudes.length)
    // }
    // for (let year = 0; year < allT[0].length; year++) {
    //     for (let lat = 0; lat < window.activeLatitudes.length; lat++) {
    //         allTbyLat[lat][year] = allT[year][lat];
    //     }
    // }
    let DATASETS_allT = [];

    for (let dim = 0; dim < allTbyLat.length; dim++) {
        DATASETS_allT.push({
            label: 'latitude: ' + Math.round(window.xLatitudes[window.activeLatitudes[dim]] * 100) / 100 + '°N',
            data: allTbyLat[dim],
            fill: false,
            borderColor: randColors[dim],
            pointRadius: 0
        })
    }
    if (window.activeLatitudes.length < 7) {
        tebm_chart_all.options.plugins.legend.display = true;
    } else {
        tebm_chart_all.options.plugins.legend.display = false;
    }
    tebm_chart_all.data.datasets = DATASETS_allT;
    tebm_chart_all.update();
}

/* DEFAULT PLOT  */
window.plot_default_tebm_chart = function plot_default_tebm_chart() {
    /* CHART WITH ONE LINE */
    document.getElementById('tempEBMchart').remove();
    document.getElementById('graph-container').innerHTML = '<canvas id="tempEBMchart"></canvas>';
    const DATA_COUNT = default_TEBM_result['x'].length;
    const LABELS = default_TEBM_result['x'].map(function (entry) {
        return Math.asin(entry) * (180 / Math.PI);
    });
    let DATASETS = [{
        label: 'Default Temperature',
        data: default_TEBM_result['T'],
        fill: false,
        borderColor: 'rgb(255,0,0)',
        pointRadius: 0
    }];

    let ctx = document.getElementById('tempEBMchart');
    let tebm_chart = new Chart(ctx, {
        type: 'line',
        data: {
            labels: LABELS,
            datasets: DATASETS,
        },
        options: {
            responsive: true,
            plugins: {
                title: {
                    display: true,
                    text: 'Temperatures at different latitudes (final year)',
                    font: {
                        Family: 'Helvetica',
                        size: 18
                    }
                },
                legend: {
                    position: 'top',
                }
            },

            scales: {
                x: {
                    display: true,
                    title: {
                        display: true,
                        text: 'Latitude',
                        font: {
                            family: 'Helvetica',
                            size: 16,
                        }
                    },
                    ticks: {
                        callback: function (value, index, values) {
                            return Math.round(LABELS[index]);
                        }
                    }
                },
                y: {
                    display: true,
                    title: {
                        display: true,
                        text: 'Temperature in °C',
                        font: {
                            family: 'Helvetica',
                            size: 16
                        }
                    }
                }
            },
            animations: {
                radius: {
                    duration: 400,
                    easing: 'linear',
                    loop: (ctx) => ctx.activate
                }
            },
            hoverRadius: 12,
            hoverBackgroundColor: 'yellow',
            interaction: {
                mode: 'nearest',
                intersect: false,
                axis: 'x'
            }

        }
    });

    window.tebm_chart = tebm_chart;

    /* CHART WITH ALL LINES */
    document.getElementById('tempEBMchart_all').remove();
    document.getElementById('tebm_all_graphs_container').innerHTML = '<canvas id="tempEBMchart_all" ></canvas>';

    let DATASETS_allT = [];
    let randColors = poolColors(default_TEBM_result['allT'].length)

    let allTbyLat = new Array(default_TEBM_result['allT'][0].length)
    for (let i = 0; i < allTbyLat.length; i++) {
        allTbyLat[i] = new Array(default_TEBM_result['allT'].length)
    }
    for (let year = 0; year < default_TEBM_result['allT'][0].length; year++) {
        for (let lat = 0; lat < default_TEBM_result['allT'].length; lat++) {
            allTbyLat[lat][year] = default_TEBM_result['allT'][year][lat];
        }
    }

    for (let dim = 0; dim < allTbyLat.length; dim++) {
        DATASETS_allT.push({
            data: allTbyLat[dim],
            fill: false,
            borderColor: randColors[dim],
            pointRadius: 0
        })
    }

    let ctx2 = document.getElementById('tempEBMchart_all');
    let tebm_chart_all = new Chart(ctx2, {
        type: 'line',
        data: {
            labels: default_TEBM_result['t'],
            datasets: DATASETS_allT,
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            plugins: {
                title: {
                    display: true,
                    text: 'Temperatures over the years at different latitudes (0°N to 90°N)',
                    font: {
                        Family: 'Helvetica',
                        size: 18
                    }
                },
                legend: {
                    display: false,
                }
            },
            scales: {
                x: {
                    display: true,
                    title: {
                        display: true,
                        text: 'Year',
                        font: {
                            family: 'Helvetica',
                            size: 16,
                        }
                    },
                    ticks: {
                        callback: function (value, index, values) {
                            return Math.round(default_TEBM_result['t'][index]);
                        }
                    }
                },
                y: {
                    display: true,
                    title: {
                        display: true,
                        text: 'Temperature °C',
                        font: {
                            family: 'Helvetica',
                            size: 16
                        }
                    }
                }
            },
            animations: {
                radius: {
                    duration: 400,
                    easing: 'linear',
                    loop: (ctx) => ctx.activate
                }
            }
        }
    });

    window.tebm_chart_all = tebm_chart_all;
    window.activeLatitudes = new Array(window.xLatitudes.length);
    for (let i = 0; i < window.xLatitudes.length; i++) {
        window.activeLatitudes[i] = parseInt(i);
    }
}

/* FUNCTION TO PLOT SELECTED LATITUDES IN TEMPERATURE / YEAR PLOT */
window.updateTYPlot = function updateTYPlot(lat_idx) {
    if (lat_idx != "all") {
        lat_idx = parseInt(lat_idx);
    }
    // show all latitudes
    if (lat_idx == "all") {
        window.activeLatitudes = [];
        for (let i = 0; i < window.last_TEBM_res['allT'].length; i++) {
            window.activeLatitudes.push(i);
        }
        document.getElementById('xlatitudes-container').selectedIndex = 0;
    } else if (window.add_this_latitude && window.activeLatitudes.includes(lat_idx)) {
        return;
    } else { // delete all and show just the one
        window.activeLatitudes = [];
        window.activeLatitudes.push(lat_idx);
    }
    updateTEBM_charts();
}

/* LATITUDES IN T / yr Plot */
window.xLatitudes = []; // real latitudes 0 to 82
window.activeLatitudes = []; // as Index of window.xLatitudes

let default_TEBM_result = scaleTEBMResults(calcTEBM()); //calcTEBM(0.6, 193, 2.1, 9.8, 420, 240, 0.7, 0.1, 0.4, 0)
window.last_TEBM_res = default_TEBM_result;
for (let lat = 0; lat < default_TEBM_result['x'].length / 2; lat++) {
    window.xLatitudes.unshift(-(Math.asin(default_TEBM_result['x'][lat]) * (180 / Math.PI)));
    window.activeLatitudes.push(lat);
}

/* create "ALL" option in latitude select*/
let option = document.createElement("option");
option.text = "all";
option.value = "all";
document.getElementById('xlatitudes-container').appendChild(option);
/* create other options in selection */
for (let i = 0; i < window.xLatitudes.length; i++) {
    let option = document.createElement("option");
    option.text = Math.round(window.xLatitudes[i]);
    option.value = i; //i as index | window.xLatitudes[i];
    document.getElementById('xlatitudes-container').appendChild(option);
}

// window.updateTYPlot = updateTYPlot;

/* global stuff */
window.onload = window.plot_default_tebm_chart();
window.added_tebm_graphs = 0;
// window.plot_default_tebm_chart = plot_default_tebm_chart;
// window.doTEBM = doTEBM;
window.default_TEBM_input = {
    D: 0.6,
    A: 193,
    B: 2.1,
    cw: 9.8,
    S0: 420,
    S2: 240,
    a0: 0.7,
    a2: 0.1,
    ai: 0.4,
    F: 0,
    gamma: 1
}

/*###################################################################################################################################################*/
/*###################################################################################################################################################*/
/*###################################################################################################################################################*/
/*###################################################################################################################################################*/
/*###################################################################################################################################################*/
/*###################################################################################################################################################*/
/*###################################################################################################################################################*/
/*###################################################################################################################################################*/
/*###################################################################################################################################################*/
/*###################################################################################################################################################*/
/*###################################################################################################################################################*/
/*###################################################################################################################################################*/
/*####################################################### COMPLEX ####################################################################*/

// D = 0.6 # diffusivity for heat transport (W m^-2 K^-1) | typical range: 0.44 - 0.66
// S1 = 338; # insolation seasonal dependence (W m^-2)
// A = 193 # OLR when T = T_m (W m^-2)
// B = 2.1 # OLR temperature dependence (W m^-2 K^-1)
// cw = 9.8 # ocean mixed layer heat capacity (W yr m^-2 K^-1) | Wärmekapazität der gemischten Ozeanschicht
// S0 = 420 # insolation at equator (W m^-2)
// S2 = 240 # insolation spatial dependence (W m^-2)
// a0 = 0.7 # ice-free coalbedo at equator
// a2 = 0.1 # ice-free coalbedo spatial dependence
// ai = 0.4 # co-albedo where there is sea ice
// Fb = 4; # heat flux from ocean below (W m^-2)
// k = 2; # sea ice thermal conductivity (W m^-2 K^-1)
// Lf = 9.5; # sea ice latent heat of fusion (W yr m^-3) | latente Schmelzwärme vom Meereis
// cg = 0.01*cw; # ghost layer heat capacity(W yr m^-2 K^-1)
// tau = 1e-5; # ghost layer coupling timescale (yr)
// winter = 26 #time of coldest <T>
// summer = 76 #time of warmest <T>

window.default_complex_ebm_input = {
    D: 0.6,
    S1: 338,
    A: 193,
    B: 2.1,
    cw: 9.8,
    S0: 420,
    S2: 240,
    a0: 0.7,
    a2: 0.1,
    ai: 0.4,
    Fb: 4,
    k: 2,
    Lf: 9.5,
    cg: 0.01,
    tau: 1e-5,
    winter: 26,
    summer: 76,
    years: 30
}

function calculate_complex_ebm(D = 0.6, S1 = 338, A = 193, B = 2.1, cw = 9.8, S0 = 420, S2 = 240,
    a0 = 0.7, a2 = 0.1, ai = 0.4, Fb = 4, k = 2, Lf = 9.5, cg = 0.01,
    tau = 1e-5, winter = 26, summer = 76, years = 30) {

    /* progressbar */
    // let progressbar = document.getElementById('progress_bar');
    // progressbar.style.width = "0%";
    // progressbar.style.display = "";
    /*  */
    D = parseFloat(D);
    S1 = parseInt(S1);
    A = parseInt(A);
    B = parseFloat(B);
    cw = parseFloat(cw);
    S0 = parseInt(S0);
    S2 = parseInt(S2);
    a0 = parseFloat(a0);
    a2 = parseFloat(a2);
    ai = parseFloat(ai);
    Fb = parseInt(Fb);
    k = parseInt(k);
    Lf = parseFloat(Lf);
    cg = parseFloat(cg);
    tau = parseFloat(tau);
    winter = parseInt(winter);
    summer = parseInt(summer);
    years = parseInt(years);

    cg *= cw;
    // #The default run in WE15, Fig 2 uses the time-stepping parameters: -------
    // n=400; % # of evenly spaced latitudinal gridboxes (equator to pole)
    // nt=1e3; % # of timesteps per year (approx lower limit of stability)
    // dur=200; % # of years for the whole run
    // #For a quicker computation, use the parameters: --------------------------
    let n = 100;
    let nt = 1e3;
    let dur = years;
    let dt = 1 / nt;
    // Spatial Grid -------------------------------------------------------------
    let dx = 1.0 / n // grid box width
    let x = range(dx / 2, 1 + dx / 2, dx); // native grid
    let xb = range(dx, 1, dx);

    // ##Diffusion Operator (WE15, Appendix A) -----------------------------------
    //in python: lam=D/dx**2*(1-xb**2)
    let lam = xb.map(function (entry, index) {
        return D / Math.pow(dx, 2) * (1 - Math.pow(entry, 2));
    })

    let L1 = [];
    L1.push(0);
    for (let i = 0; i < lam.length; i++) {
        L1.push(-lam[i]);
    }

    let L2 = [];
    for (let i = 0; i < lam.length; i++) {
        L2.push(-lam[i])
    }
    L2.push(0);

    let L3 = [];
    for (let i = 0; i < L1.length; i++) {
        L3.push(-L1[i] - L2[i]);
    }

    let L3_diag = [];
    for (let i = 0; i < L3.length; i++) {
        L3_diag.push(new Array(L3.length));
    }

    for (let row = 0; row < L3_diag.length; row++) {
        for (let column = 0; column < L3_diag.length; column++) {
            if (row === column) {
                L3_diag[row][column] = L3[row];
            } else {
                L3_diag[row][column] = 1;
                L3_diag[row][column] = 0;
            }
        }
    }

    let L2_diag = [];
    for (let i = 0; i < L2.length; i++) {
        L2_diag.push(new Array(L2.length))
    }

    for (let row = 0; row < L2_diag.length; row++) {
        for (let column = 0; column < L2_diag.length; column++) {
            if (row + 1 === column) {
                L2_diag[row][column] = L2[row];
            } else {
                L2_diag[row][column] = 1;
                L2_diag[row][column] = 0;
            }
        }
    }

    let L1_diag = [];
    for (let i = 0; i < L1.length; i++) {
        L1_diag.push(new Array(L1.length))
    }

    for (let row = 0; row < L1_diag.length; row++) {
        for (let column = 0; column < L1_diag.length; column++) {
            if (row === column + 1) {
                L1_diag[row][column] = L1[row];
            } else {
                L1_diag[row][column] = 1;
                L1_diag[row][column] = 0;
            }
        }
    }

    let diffop = [];
    for (let i = 0; i < L1.length; i++) {
        diffop.push(new Array(L1.length))
    }

    for (let row = 0; row < L1.length; row++) {
        for (let column = 0; column < L1.length; column++) {
            diffop[row][column] = -L3_diag[row][column] - L2_diag[row][column] - L1_diag[row][column];
        }
    }

    // ##Definitions for implicit scheme on Tg
    let cg_tau = cg / tau;
    let dt_tau = dt / tau;
    let dc = dt_tau * cg_tau;

    // kappa = (1+dt_tau)*np.identity(n)-dt*diffop/cg;
    let kappa = new Array(diffop.length);
    for (let i = 0; i < kappa.length; i++) {
        kappa[i] = new Array(diffop[i].length);
    }
    let I = [];
    for (let i = 0; i < n; i++) {
        I[i] = new Array(n);
    }
    for (let row = 0; row < I.length; row++) {
        for (let column = 0; column < I[row].length; column++) {
            if (row == column) {
                I[row][column] = 1;
            } else {
                I[row][column] = 1;
                I[row][column] = 0;
            }
        }
    }

    for (let row = 0; row < kappa.length; row++) {
        for (let column = 0; column < kappa[row].length; column++) {
            kappa[row][column] = (1 + dt_tau) * I[row][column] - dt * diffop[row][column] / cg;
        }
    }

    // ##Seasonal forcing (WE15 eq.3)
    let ty = range(dt / 2, 1 + dt / 2, dt)

    // S = (np.tile(S0-S2*x**2,[int(nt),1])- np.tile(S1*np.cos(2*np.pi*ty),[n,1]).T*np.tile(x,[int(nt),1]));
    // S = S___1 - S___2 * S___3
    // in python: S0-S2*x**2
    let S0S2x2 = new Array(x.length);
    for (let i = 0; i < S0S2x2.length; i++) {
        S0S2x2[i] = S0 - S2 * Math.pow(x[i], 2);
    }
    // np.tile(S0-S2*x**2,[int(nt),1])
    let S___1 = new Array(parseInt(nt));
    for (let i = 0; i < S___1.length; i++) {
        S___1[i] = S0S2x2;
    }

    // np.tile(S1*np.cos(2*np.pi*ty),[n,1]).T
    let S1cos2pity = new Array(ty.length);
    for (let i = 0; i < S1cos2pity.length; i++) {
        S1cos2pity[i] = S1 * Math.cos(ty[i] * Math.PI * 2);
    }
    let S___2 = new Array(n);
    for (let i = 0; i < S___2.length; i++) {
        S___2[i] = S1cos2pity;
    }
    S___2 = math.transpose(S___2);

    // np.tile(x,[int(nt),1])
    let S___3 = new Array(parseInt(nt));
    for (let i = 0; i < S___3.length; i++) {
        S___3[i] = x;
    }

    let S = new Array(S___1.length);
    for (let i = 0; i < S.length; i++) {
        S[i] = new Array(S___1[i].length);
    }

    for (let row = 0; row < S.length; row++) {
        for (let column = 0; column < S[row].length; column++) {
            S[row][column] = S___1[row][column] - S___2[row][column] * S___3[row][column];
        }
    }

    // ##Further definitions
    let M = B + cg_tau;
    let aw = new Array(x.length);
    for (let i = 0; i < aw.length; i++) {
        aw[i] = a0 - a2 * Math.pow(x[i], 2); //# open water albedo
    }
    let kLf = k * Lf;

    // #Set up output arrays, saving 100 timesteps/year

    let E100 = new Array(n);
    let T100 = new Array(n);
    for (let i = 0; i < n; i++) {
        E100[i] = new Array(dur * 100);
        T100[i] = new Array(dur * 100);
    }

    let p = -1;
    let m = -1;

    // #Initial conditions ------------------------------------------------------
    let T = new Array(x.length);
    for (let i = 0; i < T.length; i++) {
        T[i] = 7.5 + 20 * (1 - 2 * Math.pow(x[i], 2));
    }
    let Tg = new Array(T.length);
    for (let i = 0; i < Tg.length; i++) {
        Tg[i] = T[i];
    }
    let E = new Array(T.length);
    for (let i = 0; i < E.length; i++) {
        E[i] = T[i] * cw;
    }

    // #Integration (see WE15_NumericIntegration.pdf)----------------------------
    // #Loop over Years ---------------------------------------------------------
    for (let year = 0; year < dur; year++) {
        // #Loop within One Year-------------------------------------------------

        for (let i = 0; i < parseInt(nt); i++) {
            m = m + 1;
            // #store 100 timesteps per year
            if ((p + 1) * 10 == m) {
                p = p + 1;
            }

            for (let row = 0; row < E100.length; row++) {
                E100[row][p] = E[row];
                T100[row][p] = T[row];
            }
            // #forcing
            let awE = new Array(E.length);
            for (let entry = 0; entry < awE.length; entry++) {
                if (E[entry] > 0) {
                    awE[entry] = aw[entry];
                } else {
                    awE[entry] = 0;
                }
            }

            let aiE = new Array(E.length);
            for (let entry = 0; entry < aiE.length; entry++) {
                if (E[entry] < 0) {
                    aiE[entry] = ai;
                } else {
                    aiE[entry] = 0;
                }
            }


            //in python: alpha = aw*(E>0) + ai*(E<0) #WE15, eq.4
            let alpha = new Array(awE.length);
            for (let entry = 0; entry < alpha.length; entry++) {
                alpha[entry] = awE[entry] + aiE[entry];
            }

            let C = new Array(Tg.length);
            for (let entry = 0; entry < C.length; entry++) {
                C[entry] = alpha[entry] * S[i][entry] + cg_tau * Tg[entry] - A;
            }

            // #surface temperature
            // in python: T0 = C/(M-kLf/E) #WE15, eq.A3

            let T0 = new Array(C.length);
            for (let entry = 0; entry < T0.length; entry++) {
                T0[entry] = C[entry] / (M - kLf / E[entry]);
            }

            //in python: T = E/cw*(E>=0)+T0*(E<0)*(T0<0); #WE15, eq.9
            //  EcsE0+T0ET0
            let EcwE0 = new Array(E.length);
            for (let entry = 0; entry < EcwE0.length; entry++) {
                if (E[entry] >= 0) {
                    EcwE0[entry] = E[entry] / cw;
                } else {
                    EcwE0[entry] = 0;
                }
            }

            let T0ET0 = new Array(E.length);
            for (let entry = 0; entry < T0ET0.length; entry++) {
                if (E[entry] < 0 && T0[entry] < 0) {
                    T0ET0[entry] = T0[entry];
                } else {
                    T0ET0[entry] = 0;
                }
            }

            for (let entry = 0; entry < T.length; entry++) {
                T[entry] = EcwE0[entry] + T0ET0[entry];
            }

            // #Forward Euler on E
            // in python: E = E+dt*(C-M*T+Fb); #WE15, eq.A2

            for (let entry = 0; entry < E.length; entry++) {
                E[entry] = E[entry] + dt * (C[entry] - M * T[entry] + Fb);
            }

            // #Implicit Euler on Tg
            // in python: Tg = np.linalg.solve(kappa-np.diag(dc/(M-kLf/E)*(T0<0)*(E<0)), Tg+(dt_tau*(E/cw*(E>=0)+(ai*S[i,:]-A)/(M-kLf/E)*(T0<0)*(E<0))))
            // Tg = linear.solve(kappa-diag_dcMkLfET0E,Tg+(dt_tauEcwE+aiSA)/(MkLfET0E))
            // da hab ich ja bock drauf

            let dcMkLfET0E = new Array(E.length);
            for (let entry = 0; entry < dcMkLfET0E.length; entry++) {
                if (T0[entry] < 0 && E[entry] < 0) {
                    dcMkLfET0E[entry] = dc / (M - kLf / E[entry]);
                } else {
                    dcMkLfET0E[entry] = 0;
                }
            }

            let diag_dcMkLfET0E = new Array(dcMkLfET0E.length);
            for (let entry = 0; entry < diag_dcMkLfET0E.length; entry++) {
                diag_dcMkLfET0E[entry] = new Array(diag_dcMkLfET0E.length);
                for (let column = 0; column < diag_dcMkLfET0E[entry].length; column++) {
                    diag_dcMkLfET0E[entry][column] = 0;
                }
            }

            for (let entry = 0; entry < dcMkLfET0E.length; entry++) {
                diag_dcMkLfET0E[entry][entry] = dcMkLfET0E[entry];
            }

            let firstLinearSolvePart = new Array(kappa.length);
            for (let entry = 0; entry < firstLinearSolvePart.length; entry++) {
                firstLinearSolvePart[entry] = new Array(firstLinearSolvePart.length);
            }

            for (let row = 0; row < firstLinearSolvePart.length; row++) {
                for (let column = 0; column < firstLinearSolvePart[row].length; column++) {
                    firstLinearSolvePart[row][column] = kappa[row][column] - diag_dcMkLfET0E[row][column];
                }
            }

            /* second part of solver */
            // Tg+(dt_tau*(E/cw*(E>=0)+(ai*S[i,:]-A)/(M-kLf/E)*(T0<0)*(E<0)))
            // Tg+(dt_tau*(EcwE0+(aiSA)/(MkLfET0E0)))
            let EcwE0_ = new Array(E.length);
            for (let entry = 0; entry < EcwE0_.length; entry++) {
                if (E[entry] >= 0) {
                    EcwE0_[entry] = E[entry] / cw;
                } else {
                    EcwE0_[entry] = 0;
                }
            }
            let aiSA = new Array(S[i].length);
            for (let entry = 0; entry < aiSA.length; entry++) {
                aiSA[entry] = parseFloat(ai) * parseFloat(S[i][entry]) - A;
            }
            let MkLfET0E0 = new Array(E.length);
            for (let entry = 0; entry < MkLfET0E0.length; entry++) {
                if (T0[entry] < 0 && E[entry] < 0) {
                    MkLfET0E0[entry] = M - kLf / E[entry];
                } else {
                    MkLfET0E0[entry] = 0;
                }
            }

            // in python: dt_tau*(E/cw*(E>=0)+(ai*S[i,:]-A)/(M-kLf/E)*(T0<0)*(E<0))
            let secondLinearSolvePart = new Array(aiSA.length);
            for (let entry = 0; entry < secondLinearSolvePart.length; entry++) {
                secondLinearSolvePart[entry] = aiSA[entry] / MkLfET0E0[entry];
                if (secondLinearSolvePart[entry] == -Infinity || secondLinearSolvePart[entry] == Infinity) {
                    secondLinearSolvePart[entry] = 0;
                }
                secondLinearSolvePart[entry] += EcwE0_[entry];
                secondLinearSolvePart[entry] *= dt_tau;
                secondLinearSolvePart[entry] += Tg[entry];
            }

            // let absolutResult = linear.solve(firstLinearSolvePart, secondLinearSolvePart);
            let absolutResult = linear.solve(firstLinearSolvePart, secondLinearSolvePart);
            for (let entry = 0; entry < absolutResult.length; entry++) {
                Tg[entry] = absolutResult[entry];
            }
            // console.log("Tg after: " + Tg);
            if (i == parseInt(nt) - 1) {
                console.log("YEAR: " + year);
                let progress = year / dur * 100;
                console.log(progress + "%");
                let progressbar = window.document.getElementById('progress_bar');
                progressbar.style.width = (Math.round(progress) + "%").toString();
            }
        }
    }

    // #output only converged, final year
    let tfin = linespace(0, 1, 100);
    let Efin = new Array(E100.length);
    let Tfin = new Array(E100.length);
    for (let entry = 0; entry < E100.length; entry++) {
        Efin[entry] = new Array(E100.length);
        Tfin[entry] = new Array(E100.length);
    }
    // final year = dur*entriesPerYear-entriesPerYear
    //  e.g.: 30years*100timesteps = 3000timesteps 
    //  => 3000timesteps - 100timesteps=> timestep 2900 is the beginning of the final year
    for (let row = 0; row < E100.length; row++) {
        let TE_fin_column_idx = 0;
        for (let column = dur * E100.length - E100.length; column < E100[row].length; column++) {
            Efin[row][TE_fin_column_idx] = E100[row][column];
            Tfin[row][TE_fin_column_idx] = T100[row][column];
            TE_fin_column_idx++
        }
    }

    // # ------------------------------------------------------------------------
    // #WE15, Figure 2: Default Steady State Climatology ------------------------
    // # ------------------------------------------------------------------------

    // #compute seasonal ice edge
    let xi = new Array(100);
    for (let entry = 0; entry < xi.length; entry++) {
        xi[entry] = 0;
    }
    for (let j = 0; j < tfin.length; j++) {
        for (let entry = 0; entry < E.length; entry++) {
            E[entry] = Efin[entry][j]; // entry == row in Efin
        }

        if (E.some(isLessThanZero)) {
            for (let entry = 0; entry < E.length; entry++) {
                if (E[entry] < 0) {
                    xi[j] = x[entry];
                    break;
                }
            }
        } else {
            xi[j] = Math.max(x);
        }
    }
    // console.log("x: "+ x);
    // console.log("xi: " + xi);
    // console.log("tfin: " + tfin);
    // console.log("Tfin: ");
    // console.log(Tfin);
    // console.log("Efin: ");
    // console.log(Efin);
    // console.log("Lf: " + Lf);
    // console.log("winter: " + winter);
    // console.log("summer: " + summer);

    return {
        'x': x,
        'xi': xi,
        'tfin': tfin,
        'Tfin': Tfin,
        'Efin': Efin,
        'Lf': Lf,
        'winter': winter,
        'summer': summer
    };

}

// window.complex_ebm_default_result = calculate_complex_ebm();

/*
################################################################################################################################
########################################################################################################
##################################################################################
################### PLOTTING FUNCTIONS FOR COMPLEX EBM
*/



function complex_ebm_plot() {
    let x = window.complex_ebm_result['x'];
    let xi = window.complex_ebm_result['xi'];
    let tfin = window.complex_ebm_result['tfin'];
    let Tfin = window.complex_ebm_result['Tfin'];
    let Efin = window.complex_ebm_result['Efin']
    let Lf = window.complex_ebm_result['Lf']
    let winter = window.complex_ebm_result['winter']
    let summer = window.complex_ebm_result['summer']

    /* a) SEASONAL ENTHALPY PLOT FINAL YEAR*/
    let data_seas_enthalpy = [{
        z: Efin,
        x: tfin,
        y: x,
        type: 'contour',
        // colorbar: {
        //     title: '\(E(Jm^{-1})\)',
        //     titleside: 'right',
        //     titlefont: {
        //         size: 14,
        //         family: 'Helvetica'
        //     }
        // }
    }];

    let layout_seas_enthalpy = {
        title: 'a) \(E(Jm^{-1})\)',
        xaxis: {
            title: {
                text: 't (final year)',
            },
            linecolor: 'black',
            linewidth: 1,
            mirror: true
        },
        yaxis: {
            title: {
                text: 'x',
            },
            linecolor: 'black',
            linewidth: 1,
            mirror: true
        }
    };

    Plotly.newPlot('complex_ebm_seas_enthalpy_plot', data_seas_enthalpy, layout_seas_enthalpy);

    /* b) SEASONAL Temperature PLOT FINAL YEAR*/
    let data_seas_T = [{
        z: Tfin,
        x: tfin,
        y: x,
        type: 'contour',
        // colorbar: {
        //     title: 'Temperature in °C',
        //     titleside: 'right',
        //     titlefont: {
        //         size: 14,
        //         family: 'Helvetica'
        //     }
        // }
    }];

    let layou_seas_T = {
        title: 'b) Temperature (°C)',
        xaxis: {
            title: {
                text: 't (final year)',
            },
            linecolor: 'black',
            linewidth: 1,
            mirror: true
        },
        yaxis: {
            title: {
                text: 'x',
            },
            linecolor: 'black',
            linewidth: 1,
            mirror: true
        }
    };
    Plotly.newPlot('complex_ebm_seas_T_plot', data_seas_T, layou_seas_T);

    /* c) SEASONAL SEA ICE THICKNESS PLOT FINAL YEAR*/
    let hfin = new Array(Efin.length);
    for (let row = 0; row < hfin.length; row++) {
        hfin[row] = new Array(hfin.length);
        for (let column = 0; column < hfin[row].length; column++) {
            hfin[row][column] = (Efin[row][column] < 0) ? -Efin[row][column] / Lf : 0
        }
    }


    // get maximal value for contour scale
    // let maxContour = Math.max.apply(Math, hfin[0]);
    // for (let row = 0; row < hfin.length; row++) {
    //     let val = Math.max.apply(Math, hfin[row]);
    //     maxContour = (maxContour < val) ? val : maxContour;
    // }
    let data_seas_SeaIce = [{
        z: hfin,
        x: tfin,
        y: x,
        type: 'contour',
        // contours: {
        //     start: 1,
        //     end: maxContour,
        //     size: 1
        // }
        colorscale: [
            [0, 'rgb(255,255,255)'],
            [0.25, 'rgb(166,206,227)'],
            [0.45, 'rgb(178,223,138)'],
            [0.65, 'rgb(65,105,225)'],
            [0.85, 'rgb(0,0,205)'],
            [1, 'rgb(0,0,255']
        ]

        // colorbar: {
        //     title: 'h(m)',
        //     titleside: 'right',
        //     titlefont: {
        //         size: 14,
        //         family: 'Helvetica'
        //     }
        // }
    }];
    // THERE SHOULD BE AN OTHER LINE TO SHOW THE SEA ICE BORDER
    let data_seas_SeaIce_layout = {
        title: 'c) Sea ice thickness h(m)',
        xaxis: {
            title: {
                text: 't (final year)',
            },
            linecolor: 'black',
            linewidth: 1,
            mirror: true
        },
        yaxis: {
            title: {
                text: 'x',
            },
            linecolor: 'black',
            linewidth: 1,
            mirror: true
        }
    };
    Plotly.newPlot('complex_ebm_seas_SeaIce_plot', data_seas_SeaIce, data_seas_SeaIce_layout);


    // ##################################################################################################
    /* d) SURFACE TEMPERATURE FINAL YEAR */
    document.getElementById('complex_tsurf_graph').remove();
    document.getElementById('complex_tsurf_graph_container').innerHTML = '<canvas id="complex_tsurf_graph"></canvas>';

    const LABELS = x.map(function (entry) {
        return Math.asin(entry) * (180 / Math.PI);
    });

    let Tsummer = new Array(Tfin.length);
    for (let row = 0; row < Tsummer.length; row++) {
        Tsummer[row] = Tfin[row][summer];
    }
    let Twinter = new Array(Tfin.length);
    for (let row = 0; row < Twinter.length; row++) {
        Twinter[row] = Tfin[row][winter];
    }

    let zeroLine = new Array(Tfin.length);
    for (let entry = 0; entry < zeroLine; entry++) {
        zeroLine[entry] = 1;
        zeroLine[entry] = 0;
    }

    let DATASETS_T = [{
            label: 'summer',
            data: Tsummer,
            fill: false,
            borderColor: 'rgb(255,0,0)',
            pointRadius: 0
        },
        {
            label: 'winter',
            data: Twinter,
            fill: false,
            borderColor: 'rgb(0,0,205)',
            pointRadius: 0
        },
        {
            data: zeroLine,
            label: 'zeroLine',
            fill: false,
            borderWidth: 2,
            borderColor: 'black',
            pointRadius: 0,

        }
    ];

    let ctx_tsurf = document.getElementById('complex_tsurf_graph');
    let complex_tsurf_chart = new Chart(ctx_tsurf, {
        type: 'line',
        data: {
            labels: LABELS,
            datasets: DATASETS_T,
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            plugins: {
                title: {
                    display: true,
                    text: 'Surface Temperature',
                    font: {
                        Family: 'Helvetica',
                        size: 18
                    }
                },
                legend: {
                    position: 'top',
                    labels: {
                        filter: function (item, chart) {
                            // Logic to remove a particular legend item goes here
                            return !item.text.includes('zerozeroLine');
                        }
                    }
                }
            },
            scales: {
                x: {
                    display: true,
                    title: {
                        display: true,
                        text: 'Latitude',
                        font: {
                            family: 'Helvetica',
                            size: 16,
                        }
                    },
                    ticks: {
                        callback: function (value, index, values) {
                            return Math.round(LABELS[index]);
                        }
                    }
                },
                y: {
                    display: true,
                    title: {
                        display: true,
                        text: 'Temperature in °C',
                        font: {
                            family: 'Helvetica',
                            size: 16
                        }
                    }
                }
            },
            animations: {
                radius: {
                    duration: 400,
                    easing: 'linear',
                    loop: (ctx) => ctx.activate
                }
            },
            hoverRadius: 6,
            hoverBackgroundColor: 'yellow',
            interaction: {
                mode: 'nearest',
                intersect: false,
                axis: 'x'
            }

        }
    });

    /* e) ICE THICKNESS h when 0.7 < x < 1 FINAL YEAR */
    document.getElementById('complex_iceThickness_graph').remove();
    document.getElementById('complex_iceThickness_graph_container').innerHTML = '<canvas id="complex_iceThickness_graph"></canvas>';

    let IceThicknesssummer = new Array(hfin.length);
    for (let row = 0; row < Tsummer.length; row++) {
        IceThicknesssummer[row] = hfin[row][summer];
    }
    let IceThicknesswinter = new Array(hfin.length);
    for (let row = 0; row < Twinter.length; row++) {
        IceThicknesswinter[row] = hfin[row][winter];
    }

    let LABELSgreaterEq07 = [];
    for (let entry = 0; entry < x.length; entry++) {
        if (x[entry] >= 0.7) {
            LABELSgreaterEq07.push(Math.round(Math.asin(x[entry]) * (180 / Math.PI)));
        } else {
            IceThicknesssummer.shift();
            IceThicknesswinter.shift();
        }
    }
    let DATASETS_iceThickness = [{
            label: 'summer',
            data: IceThicknesssummer,
            fill: false,
            borderColor: 'rgb(255,0,0)',
            pointRadius: 0
        },
        {
            label: 'winter',
            data: IceThicknesswinter,
            fill: false,
            borderColor: 'rgb(0,0,205)',
            pointRadius: 0
        },
        // {
        //     data: zeroLineGEq07,
        //     fill: false,
        //     borderColor: 'rgb(0,0,0)',
        //     pointRadius: 0
        // }
    ];

    let ctx_iceThickness = document.getElementById('complex_iceThickness_graph');
    let complex_iceThickness_chart = new Chart(ctx_iceThickness, {
        type: 'line',
        data: {
            labels: LABELSgreaterEq07,
            datasets: DATASETS_iceThickness,
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            plugins: {
                title: {
                    display: true,
                    text: 'Ice thickness h when 0.7 < 1',
                    font: {
                        Family: 'Helvetica',
                        size: 18
                    }
                },
                legend: {
                    position: 'top',
                }
            },

            scales: {
                x: {
                    display: true,
                    title: {
                        display: true,
                        text: 'Latitude',
                        font: {
                            family: 'Helvetica',
                            size: 16,
                        }
                    }
                },
                y: {
                    display: true,
                    title: {
                        display: true,
                        text: 'h(m)',
                        font: {
                            family: 'Helvetica',
                            size: 16
                        }
                    }
                }
            },
            animations: {
                radius: {
                    duration: 400,
                    easing: 'linear',
                    loop: (ctx) => ctx.activate
                }
            },
            hoverRadius: 6,
            hoverBackgroundColor: 'yellow',
            interaction: {
                mode: 'nearest',
                intersect: false,
                axis: 'x'
            }

        }
    });

    /* f) SEASONAL CYCLE OF ICE THICKNESS AT THE POLE H_p */
    document.getElementById('complex_seasCycleIceThickness_graph').remove();
    document.getElementById('complex_seasCycleIceThickness_graph_container').innerHTML = '<canvas id="complex_seasCycleIceThickness_graph"></canvas>';

    let seas_iceThicknesAtPole = new Array(hfin.length);
    for (let entry = 0; entry < hfin.length; entry++) {
        seas_iceThicknesAtPole[entry] = hfin[hfin.length - 1][entry];
    }

    let LABEL_X_year = new Array(hfin.length);
    for (let entry = 0; entry < LABEL_X_year.length; entry++) {
        LABEL_X_year[entry] = entry / 100
    }
    let DATASETS_seas_iceThickness = [{
        label: 'ice tickness at pole',
        data: seas_iceThicknesAtPole,
        fill: false,
        borderColor: 'rgb(0,0,0)',
        pointRadius: 0
    }];

    let ctx_seas_iceThickness = document.getElementById('complex_seasCycleIceThickness_graph');
    let complex_seas_iceThickness_chart = new Chart(ctx_seas_iceThickness, {
        type: 'line',
        data: {
            labels: LABEL_X_year, //xi
            datasets: DATASETS_seas_iceThickness,
        },
        options: {
            responsive: true,
            plugins: {
                title: {
                    display: true,
                    text: 'Seasonal ice thickness at the pole',
                    font: {
                        Family: 'Helvetica',
                        size: 18
                    }
                },
                legend: {
                    display: false,
                }
            },
            scales: {
                x: {
                    display: true,
                    title: {
                        display: true,
                        text: 't (yr)',
                        font: {
                            family: 'Helvetica',
                            size: 16,
                        }
                    },
                    // ticks: {
                    //     callback: function (value, index, values) {
                    //         return Math.round(value) / 100;
                    //     }
                    // }
                },
                y: {
                    display: true,
                    title: {
                        display: true,
                        text: 'h(m)',
                        font: {
                            family: 'Helvetica',
                            size: 16
                        }
                    }
                }
            },
            animations: {
                radius: {
                    duration: 400,
                    easing: 'linear',
                    loop: (ctx) => ctx.activate
                }
            },
            hoverRadius: 6,
            hoverBackgroundColor: 'yellow',
            interaction: {
                mode: 'nearest',
                intersect: false,
                axis: 'x'
            }

        }
    });

    /* g) SEASONAL CYCLE OF THE LAT OF SEA ICE EDGE */
    document.getElementById('complex_seasCycleIceEdge_graph').remove();
    document.getElementById('complex_seasCycleIceEdge_graph_container').innerHTML = '<canvas id="complex_seasCycleIceEdge_graph"></canvas>';

    let DATASETS_seas_iceEdge = [{
        label: 'summer',
        data: xi, //tfin,
        fill: false,
        borderColor: 'rgb(0,0,255)',
        pointRadius: 0
    }];

    let ctx_seas_iceEdge = document.getElementById('complex_seasCycleIceEdge_graph');
    let complex_iceEdge_chart = new Chart(ctx_seas_iceEdge, {
        type: 'line',
        data: {
            labels: LABELS, // LABELS
            datasets: DATASETS_seas_iceEdge,
        },
        options: {
            responsive: true,
            plugins: {
                title: {
                    display: true,
                    text: 'Seasonal cycle of the ice edge on latitude',
                    font: {
                        Family: 'Helvetica',
                        size: 18
                    }
                },
                legend: {
                    display: false,
                    position: 'top',
                }
            },

            scales: {
                x: {
                    display: true,
                    title: {
                        display: true,
                        text: 't (yr)',
                        font: {
                            family: 'Helvetica',
                            size: 16,
                        }
                    },
                    ticks: {
                        callback: function (value, index, values) {
                            return Math.round(value) / 100;
                        }
                    }

                },
                y: {
                    display: true,
                    title: {
                        display: true,
                        text: 'Latitude',
                        font: {
                            family: 'Helvetica',
                            size: 16
                        }
                    }
                }
            },
            animations: {
                radius: {
                    duration: 400,
                    easing: 'linear',
                    loop: (ctx) => ctx.activate
                }
            },
            hoverRadius: 12,
            hoverBackgroundColor: 'yellow',
            interaction: {
                mode: 'nearest',
                intersect: false,
                axis: 'x'
            }

        }
    });

}


/*
################################################################################################################################
########################################################################################################
##################################################################################
################### RUN COMPLEX EBM
*/


window.doTheComplexThing = function doTheComplexThing(complex_ebm_input = window.default_complex_ebm_input) {
    window.complex_ebm_result = calculate_complex_ebm(complex_ebm_input['D'], complex_ebm_input['S1'],
        complex_ebm_input['A'], complex_ebm_input['B'], complex_ebm_input['cw'], complex_ebm_input['S0'], complex_ebm_input['S2'],
        complex_ebm_input['a0'], complex_ebm_input['a2'], complex_ebm_input['ai'], complex_ebm_input['Fb'], complex_ebm_input['k'],
        complex_ebm_input['Lf'], complex_ebm_input['cg'], complex_ebm_input['tau'], complex_ebm_input['winter'], complex_ebm_input['summer'],
        complex_ebm_input['years']);
    complex_ebm_plot();
}