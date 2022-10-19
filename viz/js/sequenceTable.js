
class SequenceTable {
  constructor(seq, data, globalApplicationState) {
    this.seq = seq;
    this.data = data.groupBy(['start','end','insert']);
    this.drawTable();
  }

  drawTable(){
    let tempData = globalApplicationState.selectedStartPeak.size === 0 ? this.data : this.data.filter((d) => globalApplicationState.selectedStartPeak.has(d.start));
    tempData = this.addSecondRowsToData(tempData);
    let rowSelection = d3.select('#predictionTableBody')
    .selectAll('tr')
    .data(tempData)
    .join('tr');

    rowSelection.on('click', (event, d) =>
    {
        if (true)
        {
            let placeholder = 0;
        }
    });

    let seqSelection = rowSelection.selectAll('td')
        .data((d) => {
          let type = this.determineErrorType(d);
          let errorType = {
            type: 'text',
            class: type,
            value: d.isFirst ? type : ''
          };

          let frequencyNum = {
            type: 'text',
            class: 'plain',
            value: d.isFirst ? d.values.length : ''
          };

          let sequenceOrigin = {
            type: 'text',
            class: 'plain',
            value: d.isFirst ? 'Original Sequence:' : 'Bin Sequence:'
          };

          let seq = this.determineAlignedSequence(d);

          let sequence = {
            type: 'multicolor_text',
            insertion_class: `alignment ${type}`,
            value: seq
          };

          let dataList = [errorType, frequencyNum, sequenceOrigin, sequence];
          return dataList;

        })
        .join('td')
        .attr('class', d => d.class);

    seqSelection.selectAll('tspan').remove();
    let alignmentSelection = seqSelection.filter(d => d.type === 'multicolor_text');
    alignmentSelection
      .append('tspan')
      .attr('class', 'alignment')
      .text(d => d.value.split('.')[0]);
    alignmentSelection
      .append('tspan')
      .attr('class', d => d.insertion_class)
      .text(d => d.value.split('.')[1]);
    alignmentSelection
      .append('tspan')
      .attr('class', 'alignment')
      .text(d => d.value.split('.')[2]);


    let otherSelection = seqSelection.filter(d => d.type === 'text');
    otherSelection.text(d => d.value);

  }

  addSecondRowsToData(data){
    let tempData = [];
    for (let row of data){
      row.isFirst = true;
      tempData.push(row);
      let secondRow = {...row}
      secondRow.isFirst = false;
      tempData.push(secondRow);
    }
    return tempData;
  }


  determineErrorType(d){
    if (d.start === d.end) {return 'insertion';}
    else if (d.insert.length !== 0) { return 'mutation';}
    else {return 'deletion'};
  }

  determineAlignedSequence(d){
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
        insert = d.isFirst ? '-'.repeat(insertLength) : d.insert;
        break;
      case 'deletion':
        insert = d.isFirst ? seq.substring(first_sub_end, second_sub_start) : '-'.repeat(insertLength);
        break;
      case 'mutation':
        insert = d.isFirst ? seq.substring(first_sub_end, second_sub_start) : d.insert;
        let possibleDashes = '';
        if ((d.end-d.start > d.insert.length && !d.isFirst) || (d.end-d.start < d.insert.length && d.isFirst)){
          possibleDashes = '-'.repeat(Math.abs(d.end-d.start-d.insert.length));
        }
        insert = insert + possibleDashes;
        break;
    }
    return first_sub + '.' + insert + '.' + second_sub;

  }
}
