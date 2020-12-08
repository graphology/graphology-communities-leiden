/**
 * Graphology Leiden Utils
 * ========================
 *
 * Miscellaneous utilities used by the Leiden algorithm.
 */
var SparseMap = require('mnemonist/sparse-map');
var createRandom = require('pandemonium/random').createRandom;

function addWeightToCommunity(map, community, weight) {
  var currentWeight = map.get(community);

  if (typeof currentWeight === 'undefined')
    currentWeight = 0;

  currentWeight += weight;

  map.set(community, currentWeight);
}

function UndirectedLeidenAddenda(index, options) {
  options = options || {};

  var rng = options.rng || Math.random;
  var randomness = 'randomness' in options ? options.randomness : 0.01;

  this.index = index;
  this.random = createRandom(rng);
  this.randomness = randomness;

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
  var neighboringCommunities = this.neighboringCommunities;

  var totalNodeWeight = 0;

  var i, w;
  var ei, el, et;

  // Initializing singletons
  for (i = start; i < stop; i++) {

    // Placing node in singleton
    this.belongings[i] = i;
    this.nonSingletonClusters[i] = 0;

    this.clusterWeights[i] = index.loops[i];
    this.externalEdgeWeightPerCluster[i] = 0; // TODO: loops or not?
    totalNodeWeight += index.loops[i] / 2; // TODO: how to count loops?

    ei = index.starts[i];
    el = index.starts[i + 1];

    for (; ei < el; ei++) {
      et = index.neighborhood[ei];

      // Only considering links inside of macro community
      if (index.belongings[et] !== currentMacroCommunity)
        continue;

      w = index.weights[et];
      totalNodeWeight += w;
      this.clusterWeights[i] += w;
      this.externalEdgeWeightPerCluster[i] += w;
    }
  }

  // Random iteration over nodes
  var s, ri, ci;
  var order = stop - start;

  var degree,
      bestCommunity,
      qualityValueIncrement,
      maxQualityValueIncrement,
      totalTransformedQualityValueIncrement,
      targetCommunity,
      targetCommunityDegree,
      targetCommunityWeights;

  ri = this.random(start, stop - 1);

  for (s = start; s < stop; s++, ri++) {
    i = start + (ri % order);

    // If node is not in a singleton anymore, we can skip it
    if (this.nonSingletonClusters[i] === 1)
      continue;

    // If connectivity constraint is not satisfied, we can skip the node
    if (
      this.externalEdgeWeightPerCluster[i] <
      (this.clusterWeights[i] * (totalNodeWeight - this.clusterWeights[i]) * this.resolution)
    )

    // Removing node from its current cluster
    this.clusterWeights[i] = 0;
    this.externalEdgeWeightPerCluster[i] = 0;

    // Finding neighboring clusters (including the current singleton one)
    neighboringCommunities.clear();
    neighboringCommunities.set(i, 0);

    degree = 0;

    ei = index.starts[i];
    el = index.starts[i + 1];

    for (; ei < el; ei++) {
      et = index.neighborhood[ei];

      // NOTE: we could index not to repeat this, but I have a feeling this
      // would not justify the spent memory
      if (index.belongings[et] !== currentMacroCommunity)
        continue;

      w = index.weights[et];

      degree += w;

      addWeightToCommunity(
        neighboringCommunities,
        this.belongings[et],
        w
      );
    }

    // Checking neighboring clusters
    bestCommunity = i;
    maxQualityValueIncrement = 0;
    totalTransformedQualityValueIncrement = 0;

    for (ci = 0; ci < neighboringCommunities.size; ci++) {
      targetCommunity = neighboringCommunities.dense[ci];
      targetCommunityDegree = neighboringCommunities.vals[ci];
      targetCommunityWeights = this.clusterWeights[targetCommunity];

      // Connectivity constraint
      if (
        this.externalEdgeWeightPerCluster[targetCommunity] >=
        (targetCommunityWeights * (totalNodeWeight - targetCommunityWeights) * this.resolution)
      ) {
        qualityValueIncrement = (
          targetCommunityDegree -
          degree * targetCommunityWeights * this.resolution
        );

        if (qualityValueIncrement > maxQualityValueIncrement) {
          bestCommunity = targetCommunity;
          maxQualityValueIncrement = qualityValueIncrement;
        }

        if (qualityValueIncrement >= 0)
          totalTransformedQualityValueIncrement += Math.exp(qualityValueIncrement / this.randomness);
      }

      this.cumulativeIncrement[targetCommunity] = totalTransformedQualityValueIncrement;
    }

    // TODO: continue here...
    // console.log(i, bestCommunity)
  }
};

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

exports.addWeightToCommunity = addWeightToCommunity;
exports.UndirectedLeidenAddenda = UndirectedLeidenAddenda;