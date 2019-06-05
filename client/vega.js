self.importScripts('https://cdn.jsdelivr.net/npm/vega-embed@4.0.0/build/vega-embed.js', 'https://cdn.jsdelivr.net/npm/vega-lite@3.2.1/build/vega-lite.js', 'https://cdn.jsdelivr.net/npm/vega@5.3.5/build/vega.js')

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

// async function addBases(chrom, n, end) {
//     const r = fetchChrom(chrom, end, end + n);
//     const body = await r.json();
//     let vega = await vegaEmbed('#vis', vlSpec);
//     var changeSet = vega
//         .changeset()
//         .insert(rs)
//     res.view.change('fasta', changeSet).run();
// }

// Embed the visualization in the container with id `vis`
async function buildVega(chrom, fr, to) {
    const genom = await fetchChrom(chrom, fr, to);
    const body = await genom.json();

    const al = await fetchAlignments(chrom, fr, to);
    const albody = await al.json();

    const spec = await fetchVegaSpecs();
    const vlSpec = await spec.json();
    vlSpec["width"] = $(window).width() - 150;
    let v = await vegaEmbed('#vis', vlSpec);
    v = v.view.insert("fasta", $.merge(body, albody));
    // v.addEventListener('mousedown', async function(event, item) {
    //     console.log('Zieh mich', event, item);
    //     const n = await fetchChrom(chrom, parseInt(to) + 1, parseInt(to) + 100);
    //     const upd = await n.json();
    //
    //     const m = await fetchAlignments(chrom, parseInt(to) + 1, parseInt(to) + 100);
    //     const alupd = await m.json();
    //
    //     v.insert('fasta', $.merge(upd, alupd));
    //     to = parseInt(to)  + 100;
    // });
    v.run();
}