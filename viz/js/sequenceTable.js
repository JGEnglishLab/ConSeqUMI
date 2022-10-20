
class SequenceTable {
  constructor(seq, data, globalApplicationState) {
    this.seq = seq;
    this.data = data.groupBy(['start','end','insert']).sort((a,b) => a.values.length > b.values.length ? -1 : 1);
    this.headerData = [
        {
            sorted: false,
            ascending: false,
            key: 'nucleotideIndex'
        },
        {
            sorted: false,
            ascending: false,
            key: 'errorType',
        },
        {
            sorted: false,
            ascending: false,
            key: 'frequency',
        },
    ]

    this.attachSortHandlers();
  }

  drawTable(){

    this.updateHeaders();
    let tempData = globalApplicationState.selectedStartPeak.size === 0 ? this.data : this.data.filter((d) => globalApplicationState.selectedStartPeak.has(d.start));
    d3.select('#sequenceTable').selectAll('tbody').remove();
    let rowSelection = d3.select('#sequenceTable')
      .selectAll('tbody')
      .data(tempData)
      .join('tbody');
    let consSeqRows = rowSelection
      .append('tr')
      .attr('id', 'consensusSeq');
    let binSeqRows = rowSelection
      .append('tr')
      .attr('id', 'binSeq');
    let lineSeparators = rowSelection
      .append('hr');

    let isFirst = true;
    this.rowToCellDataTransformForSelection(consSeqRows, isFirst);
    this.rowToCellDataTransformForSelection(binSeqRows, !isFirst);

  }

  rowToCellDataTransformForSelection(selection, isFirst){
    let returnSelection = selection.selectAll('td')
        .data((d) => {

          let nucleotideIndex = {
            type: 'text',
            class: 'plain',
            value: isFirst ? d.start : ''
          }

          let type = this.determineErrorType(d);
          let errorType = {
            type: 'text',
            class: type,
            value: isFirst ? type : ''
          };

          let frequencyNum = {
            type: 'text',
            class: 'plain',
            value: isFirst ? d.values.length : ''
          };

          let sequenceOrigin = {
            type: 'text',
            class: 'plain',
            value: isFirst ? 'Consensus Sequence:' : 'Bin Sequence:'
          };

          let seq = this.determineAlignedSequence(d, isFirst);

          let sequence = {
            type: 'multicolor_text',
            insertion_class: `alignment ${type}`,
            value: seq
          };

          let dataList = [nucleotideIndex, errorType, frequencyNum, sequenceOrigin, sequence];
          return dataList;

        })
        .join('td')
        .attr('class', d => d.class);
      returnSelection.selectAll('tspan').remove();
      let alignmentSelection = returnSelection.filter(d => d.type === 'multicolor_text');
      alignmentSelection
        .append('tspan')
        .attr('class', 'alignment')
        .text(d => d.value[0]);
      alignmentSelection
        .append('tspan')
        .attr('class', d => d.insertion_class)
        .text(d => d.value[1]);
      alignmentSelection
        .append('tspan')
        .attr('class', 'alignment')
        .text(d => d.value[2]);


        let otherSelection = returnSelection.filter(d => d.type === 'text');
        otherSelection.text(d => d.value);

  }

  setTextClassesAndValuesForSelection(selection){

  }

  determineErrorType(d){
    if (d.start === d.end) {return 'insertion';}
    else if (d.insert.length !== 0) { return 'mutation';}
    else {return 'deletion'};
  }

  determineAlignedSequence(d, isFirst){
    let type = this.determineErrorType(d);
    let seq = globalApplicationState.seq;
    let insertLength = type === 'insertion' ? d.insert.length : d.end-d.start
    let diff = 15 - insertLength;
    let first_sub_end = d.start;
    let second_sub_start = d.end;
    let first_sub_start = first_sub_end - Math.floor(diff/2);
    let second_sub_end = second_sub_start + Math.ceil(diff/2);
    let first_sub = seq.substring(first_sub_start, first_sub_end);
    let second_sub = seq.substring(second_sub_start, second_sub_end);
    let insert;

    switch (type) {
      case 'insertion':
        insert = isFirst ? '-'.repeat(insertLength) : d.insert;
        break;
      case 'deletion':
        insert = isFirst ? seq.substring(first_sub_end, second_sub_start) : '-'.repeat(insertLength);
        break;
      case 'mutation':
        insert = isFirst ? seq.substring(first_sub_end, second_sub_start) : d.insert;
        let possibleDashes = '';
        if ((d.end-d.start > d.insert.length && !isFirst) || (d.end-d.start < d.insert.length && isFirst)){
          possibleDashes = '-'.repeat(Math.abs(d.end-d.start-d.insert.length));
        }
        insert = insert + possibleDashes;
        break;
    }
    return [first_sub, insert, second_sub];

  }

  attachSortHandlers() {

       const headerRow = d3.select('#columnHeaders')
          .selectAll('.s')
          .data(this.headerData)
          .join('th')
          .on("click", (e, d) => {
             this.headerData.filter(s => s !== d).forEach(col => col.sorted = false);
             this.headerData.filter(s => s !== d).forEach(col => col.ascending = false);
             d.sorted = true;
             let newData = this.data.sort((a,b) => this.sortColumn(a,b,d));
             d.ascending = !d.ascending;
             this.data = newData;
             this.drawTable();
          });

  }

  updateHeaders() {
   d3.select('#columnHeaders')
      .selectAll('.s')
      .data(this.headerData)
      .join('th')
      .attr("class", d => d.sorted ? "s sorting" : "s sortable");

    d3.select('#columnHeaders')
       .selectAll('i')
       .data(this.headerData)
       .join('i')
       .attr("class", (d) => {
         let sortedClass = d.ascending ? 'fas fa-sort-up' : 'fas fa-sort-down';
         return d.sorted ? sortedClass : "fas no-display";
       });


}

  sortColumn(a,b,d){
    let aVal;
    let bVal;
    switch(d.key){
      case 'nucleotideIndex':
        aVal = a.start;
        bVal = b.start;
        break;
      case 'errorType':
        aVal = this.determineErrorType(a);
        bVal = this.determineErrorType(b);
        break;
      case 'frequency':
        aVal = a.values.length;
        bVal = b.values.length;
        break;
    }
      if (d.ascending){
        return aVal > bVal ? -1 : 1
      } else {
        return aVal < bVal ? -1 : 1
     }
    }


}
