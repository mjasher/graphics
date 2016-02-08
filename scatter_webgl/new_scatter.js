// select svg canvas
function scatter(el, data, xcol, ycol){
  var m = [30, 10, 10, 10],       // margins
      w = 500,                    // width
      h = 500,                    // height
      dimensions = [],            // quantitative dimensions
      // xcol = 0,                   // active x column
      // ycol = 1,                   // active y column
      last = [],                  // last [x,y,color] pairs
      xscale = d3.scale.linear(), // x scale
      yscale = d3.scale.linear(); // yscale 

  // color scale
  var color = d3.scale.category20();

  // adjust canvas size
  var canvas = el
    .attr("width", w + "px")
    .attr("height", h + "px");

  // rendering context
  WebGL2D.enable(canvas.node())
  var ctx = canvas.node().getContext('webgl-2d');
  ctx.strokeStyle = "rgba(0,0,0,0.8)";
  ctx.lineWidth = "1.5";

  // load data from csv file

  // get columns of csv, mark excluded columns
  var columns = d3.keys( data[0] );


  // extents for each dimension
  var extents = _(dimensions)
    .map(function(col) {
      return [0, d3.max(data, function(d) { return parseFloat(d[col]) })]
    });

  // create scales
  // xscale.domain(extents[xcol]).range([m[3], w - m[1]]),
  // yscale.domain(extents[ycol]).range([h - m[2], m[0]]);

  xscale.domain(d3.extent(data, function(d){ return d[xcol]; })).range([m[3], w - m[1]]),
  yscale.domain(d3.extent(data, function(d){ return d[ycol]; })).range([h - m[2], m[0]]);


  // from data point, return [x,y,color]
  function position(d) {
    var x = xscale(d[xcol]);
    var y = yscale(d[ycol]);
    return [x, y , color(d.group)];
  };

  // render circle [x,y,color]
  function circle(pos) {
    ctx.fillStyle = pos[2];
    ctx.fillRect(pos[0], pos[1], 2, 2);
  };

  // clear canvas
  ctx.clearRect(0,0,w,h);
  // render initial data points
  data.map(position).forEach(circle);
  ctx.clearRect(0,0,w,h);

}

d3.csv('nutrients.csv', function(data) {
    var cols = ["protein (g)", "calcium (g)", "sodium (g)"];

    var charts = [];
    for (var i = 0; i < cols.length; i++) {
      var row = d3.select("body").append('div').attr('class', 'row');
      charts.push([]);
      for (var j = 0; j < cols.length; j++) {
        var cell = row.append('div').attr('class', 'cell');
        cell.append("canvas");
        if (i == j) cell.append('text').text(cols[i]);
        charts[i].push(cell);
        // charts[i].push(row.append("canvas"));
      }
    }

    for (var i = 0; i < cols.length; i++) {
      for (var j = 0; j < cols.length; j++) {
        scatter(charts[i][j].select('canvas'), data, cols[i], cols[j]);
      };
    };
    console.log(Object.keys(data[0]))

});


  // // change x axis
  // function xaxis(i) {
  //   xcol = i;
  //   xscale.domain(extents[i]);
  //   draw();
  // };

  // // change y axis
  // function yaxis(i) {
  //   ycol = i;
  //   yscale.domain(extents[i]);
  //   draw();
  // };


  // create dropdowns to change axes
  // d3.select("#xaxis")
  //   .selectAll("option")
  //   .data(dimensions)
  //   .enter().append("option")
  //     .attr("value", function(d,i) { return i; })
  //     .text(function(d) { return d; })
  //     .each(function(d,i) {
  //       if (i == xcol) d3.select(this).attr("selected", "yes");
  //     });

  // d3.select("#xaxis")
  //     .on("change", function() { 
  //       xaxis(this.selectedIndex) 
  //     });

  // d3.select("#yaxis")
  //   .selectAll("option")
  //   .data(dimensions)
  //   .enter().append("option")
  //     .attr("value", function(d,i) { return i; })
  //     .text(function(d) { return d; })
  //     .each(function(d,i) {
  //       if (i == ycol) d3.select(this).attr("selected", "yes");
  //     });

  // d3.select("#yaxis")
  //     .on("change", function() { 
  //       yaxis(this.selectedIndex) 
  //     });

  // window.data = data;

