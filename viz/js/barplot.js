/** Class implementing the table. */
class Barplot {
    /**
     * Creates a Table Object
     */
  constructor() {
    this.initializeBarPlot();
    this.drawBarPlot();
  }

  initializeBarPlot(){
    d3.select('#Barchart-div')
      .append('svg')
      .attr('id', 'Barchart-svg')
      .append('g')
      .attr('id', 'Barchart-x-axis');
    d3.select('#Barchart-svg')
      .append('g')
      .attr('id', 'Barchart-y-axis');
    d3.select('#Barchart-svg')
      .append('g')
      .attr('id', 'BarChart')
      .attr('class', 'bar-chart');


  }

  drawBarPlot(){

    const binNum = d3.select('#bin_num_select').property('value');
    const data = d3.group(globalApplicationState.data.filter(d => d.binNum === parseInt(binNum)), d => d.start);

    const yScale = d3.scaleLinear()
      .domain([0, d3.max([...data.values()], n => n.length)])
      .range([CHART_HEIGHT - MARGIN.bottom - MARGIN.top, 0])
      .nice();

    const xScale = d3.scaleLinear()
      .domain([0, globalApplicationState.seq.length])
      .range([MARGIN.left, CHART_WIDTH-MARGIN.right])
      .nice();

    d3.select('#Barchart-y-axis')
      .call(d3.axisLeft(yScale))
      .attr('transform', `translate(${MARGIN.left}, ${MARGIN.top})`);

    d3.select('#Barchart-x-axis')
      .attr('transform', `translate(0,${CHART_HEIGHT - MARGIN.bottom})`)
      .call(d3.axisBottom(xScale));

    d3.select('#BarChart')
      .selectAll('rect')
      .data(data)
      .join(
        enter => enter
          .append('rect')
          .attr('width', 2)
          .attr('x', d => xScale(d[0]))
          .attr('y', d => yScale(d[1].length) + MARGIN.top)
          .attr('height', d=> yScale(0) - yScale(d[1].length))
          .attr('opacity', 0)
          .attr('id', d => d.start)
          .transition()
          .duration(ANIMATION_DURATION)
          .delay(ANIMATION_DURATION)
          .attr('height', d=> yScale(0) - yScale(d[1].length))
          .attr('opacity', 1),
        update => update
          .transition()
          .duration(ANIMATION_DURATION)
          .attr('width', 2)
          .attr('x', d => xScale(d[0]))
          .attr('y', d => yScale(d[1].length) + MARGIN.top)
          .attr('height', d=> yScale(0) - yScale(d[1].length))
          .attr('id', d => d.start),
        exit => exit
          .transition()
          .duration(ANIMATION_DURATION)
          .attr('width', 0)
          .attr('height', 0)
          .remove()
      )
      .on('mouseover', (e) => d3.select(e.target).classed('hovered', true))
      .on('mouseout', (e) => d3.select(e.target).classed('hovered', false))
      .on('click', (e,d) => {
        globalApplicationState.selectedStartPeak.has(d[0]) ? globalApplicationState.selectedStartPeak.delete(d[0]) : globalApplicationState.selectedStartPeak.add(d[0])
        let isSelected = globalApplicationState.selectedStartPeak.has(d[0]) ? true : false;
        d3.select(e.target).classed('hovered', false);
        this.updateSelectedRects();
        globalApplicationState.sequenceTable.setClassData();
        globalApplicationState.sequenceTable.drawTable();
      });
      this.updateSelectedRects();
  }

  updateSelectedRects(){
    d3.select('#BarChart')
      .selectAll('rect')
      .classed('selected', (d) => globalApplicationState.selectedStartPeak.has(d[0]) ? true : false);
  }

}
