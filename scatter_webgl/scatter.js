// select svg canvas
var m = [30, 10, 10, 10],       // margins
    w = 500,                    // width
    h = 400,                    // height
    dimensions = [],            // quantitative dimensions
    xcol = 0,                   // active x column
    ycol = 1,                   // active y column
    last = [],                  // last [x,y,color] pairs
    transition_count = 0,       // used to cancel old transitions
    xscale = d3.scale.linear(), // x scale
    yscale = d3.scale.linear(); // yscale 

// color scale
var color = {
  "Baby Foods"                        : '#555555',
  "Baked Products"                    : '#7f7f7f',
  "Beverages"                         : '#c49c94',
  "Breakfast Cereals"                 : '#9467bd',
  "Cereal Grains and Pasta"           : '#bcbd22',
  "Dairy and Egg Products"            : '#ff7f0e',
  "Ethnic Foods"                      : '#e7ba52',
  "Fast Foods"                        : '#dbdb8d',
  "Fats and Oils"                     : '#ffbb78',
  "Finfish and Shellfish Products"    : '#e377c2',
  "Fruits and Fruit Juices"           : '#c5b0d5',
  "Legumes and Legume Products"       : '#f7b6d2',
  "Meals, Entrees, and Sidedishes"    : '#17becf',
  "Nut and Seed Products"             : '#8c564b',
  "Pork Products"                     : "#00ee99",
  "Poultry Products"                  : '#d62728',
  "Restaurant Foods"                  : '#1f77b4',
  "Sausages and Luncheon Meats"       : '#ff9896',
  "Snacks"                            : '#9edae5',
  "Soups, Sauces, and Gravies"        : '#98df8a',
  "Spices and Herbs"                  : '#aec7e8',
  "Sweets"                            : '#c7c7c7',
  "Vegetables and Vegetable Products" : '#2ca02c'
};

// adjust canvas size
canvas = d3.select("#chart")
  .attr("width", w + "px")
  .attr("height", h + "px");

// rendering context
WebGL2D.enable(canvas[0][0])
ctx = canvas[0][0].getContext('webgl-2d');
ctx.strokeStyle = "rgba(0,0,0,0.8)";
ctx.lineWidth = "1.5";

// load data from csv file
d3.csv('nutrients.csv', function(data) {

  // get columns of csv, mark excluded columns
  var columns = d3.keys( data[0] ),
      excluded = ['name', 'group', 'id'];

  // get quantitative dimensions
  dimensions = _(columns)
    .difference(excluded);

  // extents for each dimension
  var extents = _(dimensions)
    .map(function(col) {
      return [0, d3.max(data, function(d) { return parseFloat(d[col]) })]
    });

  // create scales
  xscale.domain(extents[xcol]).range([m[3], w-m[1]]),
  yscale.domain(extents[ycol]).range([h-m[2], m[0]]);

  // render initial data points
  last = data.map(position)
  last.forEach(circle);

  // change x axis
  function xaxis(i) {
    xcol = i;
    xscale.domain(extents[i]);
    transition(++transition_count);
  };

  // change y axis
  function yaxis(i) {
    ycol = i;
    yscale.domain(extents[i]);
    transition(++transition_count);
  };

  // create dropdowns to change axes
  d3.select("#xaxis")
    .selectAll("option")
    .data(dimensions)
    .enter().append("option")
      .attr("value", function(d,i) { return i; })
      .text(function(d) { return d; })
      .each(function(d,i) {
        if (i == xcol) d3.select(this).attr("selected", "yes");
      });

  d3.select("#xaxis")
      .on("change", function() { xaxis(this.selectedIndex) });

  d3.select("#yaxis")
    .selectAll("option")
    .data(dimensions)
    .enter().append("option")
      .attr("value", function(d,i) { return i; })
      .text(function(d) { return d; })
      .each(function(d,i) {
        if (i == ycol) d3.select(this).attr("selected", "yes");
      });

  d3.select("#yaxis")
      .on("change", function() { yaxis(this.selectedIndex) });

  window.data = data;
});

function transition(count) {
  // next positions
  var next = data.map(position);

  var transition = d3.interpolate(last, next);

  // run transition
  d3.timer(function(t) {
    // abort old transition
    if (count < transition_count) return true;

    clear();
    if (t > 1000) {
      last = next;
      last.forEach(circle);
      return true
    };
    last = transition(t/1000);
    last.forEach(circle);
  });
};

// clear canvas
function clear() {
  ctx.clearRect(0,0,w,h);
};

// from data point, return [x,y,color]
function position(d) {
  var x = xscale(d[dimensions[xcol]]);
  var y = yscale(d[dimensions[ycol]]);
  return [x,y,color[d.group]];
};

// render circle [x,y,color]
function circle(pos) {
  ctx.fillStyle = pos[2];
  ctx.fillRect(pos[0],pos[1],2, 2);
  /*   PREVIOUS CODE
  ctx.fillStyle = pos[2];
  ctx.beginPath();
  ctx.arc(pos[0],pos[1],2,0,2*Math.PI);
  ctx.stroke();
  ctx.fill();
  */
};