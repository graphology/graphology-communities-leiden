/**
 * Graphology Leiden Utils
 * ========================
 *
 * Miscellaneous utilities used by the Leiden algorithm.
 */
var SparseMap = require('mnemonist/sparse-map');
var createRandomIndex = require('pandemonium/random-index').createRandomIndex;

function UndirectedLeidenAddenda(index, options) {
  options = options || {};

  var rng = options.rng || Math.random;

  this.index = index;
  this.randomIndex = createRandomIndex(rng);

  var NodesPointerArray = index.counts.constructor;
  var WeightsArray = index.weights.constructor;

  var order = index.C;
  this.resolution = index.resolution;

  // Used to group nodes by communities
  this.B = 0;
  this.communitiesOffsets = new NodesPointerArray(order);
  this.nodesSortedByCommunities = new NodesPointerArray(order);
  this.communitiesBounds = new NodesPointerArray(order + 1);

  // Used to merge nodes subsets
  this.clusterWeights = new WeightsArray(order);
  this.nonSingletonClusters = new Uint8Array(order);
  this.externalEdgeWeightPerCluster = new WeightsArray(order);
  this.belongings = new NodesPointerArray(order);
  this.neighboringCommunities = new SparseMap(order);
  this.cumulativeIncrement = new Float64Array(order);
}

UndirectedLeidenAddenda.prototype.groupByCommunities = function() {
  var index = this.index;

  var n, i, c, b, o;

  n = 0;
  o = 0;

  for (i = 0; i < index.C; i++) {
    c = index.counts[i];

    if (c !== 0) {
      this.communitiesBounds[o++] = n;
      n += c;
      this.communitiesOffsets[i] = n;
    }
  }

  this.communitiesBounds[o] = n;

  o = 0;

  for (i = 0; i < index.C; i++) {
    b = index.belongings[i];
    o = --this.communitiesOffsets[b];
    this.nodesSortedByCommunities[o] = i;
  }

  this.B = index.C - index.U;
};

UndirectedLeidenAddenda.prototype.communities = function() {
  var communities = new Array(this.B);

  var i, j, community, start, stop;

  for (i = 0; i < this.B; i++) {
    start = this.communitiesBounds[i];
    stop = this.communitiesBounds[i + 1];
    community = []; // NOTE: size can be given

    for (j = start; j < stop; j++) {
      community.push(j);
    }

    communities[i] = community;
  }

  return communities;
};

UndirectedLeidenAddenda.prototype.mergeNodesSubset = function(start, stop) {
  var index = this.index;

  var currentMacroCommunity = index.belongings[this.nodesSortedByCommunities[start]];

  console.log(start, stop, currentMacroCommunity);
};

// var currentCommunity = index.belongings[nodes[start]];

// // Initializing counters
// for (j = 0, i = start; i < stop; i++, j++) {
//   n = nodes[i];
//   nodesOrder[j] = n;
//   belongings[j] = j;
//   ei = index.starts[n];
//   el = index.starts[n + 1];

//   clusterWeights[j] += index.loops[n];

//   for(; ei < el; ei++) {
//     et = index.neighborhood[ei];

//     // Only considering links from the same community
//     // TODO: need to recreate a sub neighorhood here...
//     if (index.belongings[et] !== currentCommunity)
//       continue;

//     w = index.weights[et];
//     totalNodeWeight += w;
//     clusterWeights[j] += w;
//     externalEdgeWeightPerCluster[j] += w;
//   }

UndirectedLeidenAddenda.prototype.refinePartition = function() {
  this.groupByCommunities();

  var i, start, stop;

  var bounds = this.communitiesBounds;

  for (i = 0; i < this.B; i++) {
    start = bounds[i];
    stop = bounds[i + 1];

    this.mergeNodesSubset(start, stop);
  }
};

exports.UndirectedLeidenAddenda = UndirectedLeidenAddenda;
