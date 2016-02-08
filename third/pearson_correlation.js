/**
 Slightly edited from

 *  @fileoverview Pearson correlation score algorithm.
 *  @author matt.west@kojilabs.com (Matt West)
 *  @license Copyright 2013 Matt West.
 *  Licensed under MIT (http://opensource.org/licenses/MIT).
 */

function pearson_correlation(p1, p2) {

  if (p1.length != p2.length) {
      throw "p1.length != p2.length";
  }

  var n = p1.length;

  if (n == 0) return 0;

  var sum1 = 0;
  for (var i = 0; i < n; i++) sum1 += p1[i];

  var sum2 = 0;
  for (var i = 0; i < n; i++) sum2 += p2[i];

  var sum1Sq = 0;
  for (var i = 0; i < n; i++) {
    sum1Sq += Math.pow(p1[i], 2);
  }

  var sum2Sq = 0;
  for (var i = 0; i < n; i++) {
    sum2Sq += Math.pow(p2[i], 2);
  }

  var pSum = 0;
  for (var i = 0; i < n; i++) {
    pSum += p1[i] * p2[i];
  }

  var num = pSum - (sum1 * sum2 / n);
  var den = Math.sqrt((sum1Sq - Math.pow(sum1, 2) / n) *
      (sum2Sq - Math.pow(sum2, 2) / n));

  if (den == 0) return 0;

  return num / den;
}

// console.log(pearson_correlation(new Array(20,54,54,65,45), new Array(21,54,60,78,82) ));
