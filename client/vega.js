async function fetchChrom(chrom, fr, to) {
    const rs = await fetch('/api/v1/reference/' + chrom +'/' + fr + '/' + to);
    return rs;
}

async function fetchAlignments(chrom, fr, to) {
    const rs = await fetch('/api/v1/alignment/' + chrom +'/' + fr + '/' + to);
    return rs;
}

async function fetchVegaSpecs() {
    const vlSpec = await fetch( "vlSpec.json");
    return vlSpec;
}

var lastLowerBound;
var lastUpperBound;


// Embed the visualization in the container with id `vis`
async function buildVega(chrom, fr, to) {
    lastLowerBound = fr;
    lastUpperBound = to;

    const genom = await fetchChrom(chrom, fr, to);
    const body = await genom.json();

    const al = await fetchAlignments(chrom, fr, to);
    const albody = await al.json();

    const spec = await fetchVegaSpecs();
    const vlSpec = await spec.json();
    vlSpec["width"] = $(window).width() - 150;
    let v = await vegaEmbed('#vis', vlSpec);
    v = v.view.insert("fasta", $.merge(body, albody));

    v.addEventListener('mousedown', function(event, item) {
        var myState = v.getState();
        if(myState.signals){
            if(myState.signals.grid_position){
                myState.signals.grid_position[0] = 1337;
                myState.signals.grid_position[1] = 2337;
                console.log(myState.signals.grid_position);

                v.setState(myState);
            } else {
                console.log("Fehler Null")
            }
        }
    });


   v.addEventListener('mouseup', async function (event, item) {
        const lowerBound = Math.round(v.getState().signals.grid_position[0]);
        const upperBound = Math.round(v.getState().signals.grid_position[1]);

        var upd;


        if (lastUpperBound < upperBound) {
            const n = await fetchChrom(chrom, lastUpperBound, upperBound);
            const upper_upd_ref = await n.json();

            const m = await fetchAlignments(chrom, lastUpperBound, upperBound);
            const upper_upd_al = await m.json();

            upd = $.merge(upper_upd_al, upper_upd_ref);
        } else {
            const o = await fetchChrom(chrom, lowerBound, lastLowerBound);
            const lower_upd_ref = await o.json();

            const p = await fetchAlignments(chrom, lowerBound, lastLowerBound);
            const lower_upd_al = await p.json();

            upd = $.merge(lower_upd_al, lower_upd_ref);
        }

        v.change('fasta', vega.changeset().insert(upd).remove(function (d) {
            return (d.position < lowerBound) || (d.position > upperBound);
        }))



        lastLowerBound = lowerBound;
        lastUpperBound = upperBound;
    });
    v.run();


}