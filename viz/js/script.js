
setup();
// Constants for the charts, that would be useful.
const CHART_WIDTH = 500;
const CHART_HEIGHT = 250;
const MARGIN = { left: 50, bottom: 20, top: 20, right: 20 };
const ANIMATION_DURATION = 300;

/**
 * Initial setup for the page. This should only be called once.
 */
 const globalApplicationState = {
   selectedStartPeak: new Set(),
   seq:null,
   diffData:null,
   barplot: null,
   table:null
 };

function setup () {

  let seq = 'TGCCACCCGGCTTCAACGAGTACGACTTCGTGCCCGAGAGCTTCGACCGGGACAAAACCATCGCCCTGATCATGAACAGTAGTGGCAGTACCGGATTGCCCAAGGGCGTAGCCCTACCGCACCGCACCGCTTGTGTCCGATTCAGTCATGCCCGCGACCCCATCTTCGGCAACCAGATCATCCCCGACACCGCTATCCTCAGCGTGGTGCCATTTCACCACGGCTTCGGCATGTTCACCACGCTGGGCTACTTGATCTGCGGCTTTCGGGTCGTGCTCATGTACCGCTTCGAGGAGGAGCTATTCTTGCGCAGCTTGCAAGACTATAAGATTCAATCTGCCCTGCTGGTGCCCACACTATTTAGCTTCTTCGCTAAGAGCACTCTCATCGACAAGTACGACCTAAGCAACTTGCACGAGATCGCCAGCGGCGGGGCGCCGCTCAGCAAGGAGGTAGGTGAGGCCGTGGCCAAACGCTTCCACCTACCAGGCATCCGCCAGGGCTACGGCCTGACAGAAACAACCAGCGCCATTCTGATCACCCCCGAAGGGGACGACAAGCCTGGCGCAGTAGGCAAGGTGGTGCCCTTCTTCGAGGCTAAGGTGGTGGACTTGGACACCGGCAAGACACTGGGTGTGAACCAGCGCGGCGAGCTGTGCGTCCGTGGCCCCATGATCATGAGCGGCTACGTTAACAACCCCGAGGCTACAAACGCTCTCATCGACAAGGACGGCTGGCTGCACAGCGGCGACATCGCCTACTGGGACGAGGACGAGCACTTCTTCATCGTGGACCGGCTGAAGAGCCTGATCAAATACAAGGGCTACCAGGTAGCCCCAGCCGAACTGGAGAGCATCCTGCTGCAACACCCCAACATCTTCGACGCCGGGGTCGCCGGCCTGCCCGACGACGATGCCGGCGAGCTGCCCGCCGCAGTCGTCGTGCTGGAACACGGTAAAACCATGACCGAGAAGGAGATCGTGGACTATGTGGCCAGCCAGGTTACAACCGCCAAGAAGCTGCGCGGTGGTGTTGTGTTCGTGGACGAGGTGCCTAAAGGACTGACCGGCAAGTTGGACGCCCGCAAGATCCGCGAGATTCTCATTAAGGCCAAGAAGGGCGGCAAGATCGCCGTGTAATAATCGTTCCTCTAGAGACGCGCAGGAGAAATTAATCAAGACTAGTACACTCCCCGTCGATCAGGGTGGTTACGTCAGTCACCGGTCGACTGTGCCTTCTAGTTGCCAGCCATCTGTTGTTTGCCCCTCCCCCGTGCCTTCCTTGACCCTGGAAGGTGCCACTCCCACTGTCCTTTCCTAATA';

  Promise.all([d3.csv('./data/differences.csv')]).then( data => {
    let diffData = data[0];
    diffData = diffData.map(function(d) {
      d.start = parseInt(d.start);
      d.end = parseInt(d.end);
      return d;
    })
    d3.select("#num_sequences").text(`Total Number of Sequences in Bin: ${new Set(diffData.map(d => d.seqID)).size}`);
    globalApplicationState.seq = seq;
    globalApplicationState.diffData = diffData;
    let barplot = new Barplot(seq, diffData, globalApplicationState);
    let sequenceTable = new SequenceTable(seq, diffData, globalApplicationState);
    sequenceTable.drawTable();
    globalApplicationState.barplot = barplot;
    globalApplicationState.sequenceTable = sequenceTable;
  })
}


Array.prototype.groupBy = function (props) {
   var arr = this;
   var partialResult = {};

   arr.forEach(el=>{

       var grpObj = {};

       props.forEach(prop=>{
             grpObj[prop] = el[prop]
       });

       var key = JSON.stringify(grpObj);

       if(!partialResult[key]) partialResult[key] = [];

       partialResult[key].push(el);

   });

   var finalResult = Object.keys(partialResult).map(key=>{
      var keyObj = JSON.parse(key);
      keyObj.values = partialResult[key];
      return keyObj;
   })

   return finalResult;
}
