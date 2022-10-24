
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
   data:null,
   filteredData:null,
   barplot: null,
   table:null
 };

function setup () {
  document.querySelector('#bin_num_select').addEventListener('change', changeData);

  let seq = 'TGCCACCCGGCTTCAACGAGTACGACTTCGTGCCCGAGAGCTTCGACCGGGACAAAACCATCGCCCTGATCATGAACAGTAGTGGCAGTACCGGATTGCCCAAGGGCGTAGCCCTACCGCACCGCACCGCTTGTGTCCGATTCAGTCATGCCCGCGACCCCATCTTCGGCAACCAGATCATCCCCGACACCGCTATCCTCAGCGTGGTGCCATTTCACCACGGCTTCGGCATGTTCACCACGCTGGGCTACTTGATCTGCGGCTTTCGGGTCGTGCTCATGTACCGCTTCGAGGAGGAGCTATTCTTGCGCAGCTTGCAAGACTATAAGATTCAATCTGCCCTGCTGGTGCCCACACTATTTAGCTTCTTCGCTAAGAGCACTCTCATCGACAAGTACGACCTAAGCAACTTGCACGAGATCGCCAGCGGCGGGGCGCCGCTCAGCAAGGAGGTAGGTGAGGCCGTGGCCAAACGCTTCCACCTACCAGGCATCCGCCAGGGCTACGGCCTGACAGAAACAACCAGCGCCATTCTGATCACCCCCGAAGGGGACGACAAGCCTGGCGCAGTAGGCAAGGTGGTGCCCTTCTTCGAGGCTAAGGTGGTGGACTTGGACACCGGCAAGACACTGGGTGTGAACCAGCGCGGCGAGCTGTGCGTCCGTGGCCCCATGATCATGAGCGGCTACGTTAACAACCCCGAGGCTACAAACGCTCTCATCGACAAGGACGGCTGGCTGCACAGCGGCGACATCGCCTACTGGGACGAGGACGAGCACTTCTTCATCGTGGACCGGCTGAAGAGCCTGATCAAATACAAGGGCTACCAGGTAGCCCCAGCCGAACTGGAGAGCATCCTGCTGCAACACCCCAACATCTTCGACGCCGGGGTCGCCGGCCTGCCCGACGACGATGCCGGCGAGCTGCCCGCCGCAGTCGTCGTGCTGGAACACGGTAAAACCATGACCGAGAAGGAGATCGTGGACTATGTGGCCAGCCAGGTTACAACCGCCAAGAAGCTGCGCGGTGGTGTTGTGTTCGTGGACGAGGTGCCTAAAGGACTGACCGGCAAGTTGGACGCCCGCAAGATCCGCGAGATTCTCATTAAGGCCAAGAAGGGCGGCAAGATCGCCGTGTAATAATCGTTCCTCTAGAGACGCGCAGGAGAAATTAATCAAGACTAGTACACTCCCCGTCGATCAGGGTGGTTACGTCAGTCACCGGTCGACTGTGCCTTCTAGTTGCCAGCCATCTGTTGTTTGCCCCTCCCCCGTGCCTTCCTTGACCCTGGAAGGTGCCACTCCCACTGTCCTTTCCTAATA';

  Promise.all([d3.csv('./data/differences.csv')]).then( csvData => {
    let data = csvData[0];
    data = data.map(function(d) {
      d.start = parseInt(d.start);
      d.end = parseInt(d.end);
      d.binNum = parseInt(d.binNum);
      return d;
    });
    //.filter(d => (d.binNum === 0 || d.binNum === 1));

    const binNums = [...new Set(data.map(d => d.binNum))];
    binNums.sort((a,b) => a-b);
    d3.select('#bin_num_select').selectAll('option')
      .data(binNums)
      .enter()
      .append("option")
      .text(d => d)
      .attr("value",d => d);

    globalApplicationState.seq = seq;
    globalApplicationState.data = data;
    const binNum = d3.select('#bin_num_select').property('value');
    let tempData = globalApplicationState.data.filter(d => d.binNum === parseInt(binNum))
    d3.select("#num_sequences").text(`Total Number of Sequences in Bin: ${new Set(tempData.map(d => d.seqID)).size}`);


    let sequenceTable = new SequenceTable();
    globalApplicationState.sequenceTable = sequenceTable;
    let barplot = new Barplot();
    globalApplicationState.barplot = barplot;
  })
}

function changeData(){
  const binNum = d3.select('#bin_num_select').property('value');
  let tempData = globalApplicationState.data.filter(d => d.binNum === parseInt(binNum))
  d3.select("#num_sequences").text(`Total Number of Sequences in Bin: ${new Set(tempData.map(d => d.seqID)).size}`);
  globalApplicationState.barplot.drawBarPlot();
  globalApplicationState.sequenceTable.initializeHeaderData();
  globalApplicationState.sequenceTable.setClassData();
  globalApplicationState.sequenceTable.drawTable();
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
