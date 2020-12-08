/**
 * Graphology Leiden Algorithm
 * ============================
 *
 * JavaScript implementation of the Leiden community detection
 * algorithm for graphology.
 *
 * [Articles]
 * M. E. J. Newman, « Modularity and community structure in networks »,
 * Proc. Natl. Acad. Sci. USA, vol. 103, no 23,‎ 2006, p. 8577–8582
 * https://dx.doi.org/10.1073%2Fpnas.0601602103
 *
 * Newman, M. E. J. « Community detection in networks: Modularity optimization
 * and maximum likelihood are equivalent ». Physical Review E, vol. 94, no 5,
 * novembre 2016, p. 052315. arXiv.org, doi:10.1103/PhysRevE.94.052315.
 * https://arxiv.org/pdf/1606.02319.pdf
 *
 * Blondel, Vincent D., et al. « Fast unfolding of communities in large
 * networks ». Journal of Statistical Mechanics: Theory and Experiment,
 * vol. 2008, no 10, octobre 2008, p. P10008. DOI.org (Crossref),
 * doi:10.1088/1742-5468/2008/10/P10008.
 * https://arxiv.org/pdf/0803.0476.pdf
 *
 * Nicolas Dugué, Anthony Perez. Directed Louvain: maximizing modularity in
 * directed networks. [Research Report] Université d’Orléans. 2015. hal-01231784
 * https://hal.archives-ouvertes.fr/hal-01231784
 *
 * R. Lambiotte, J.-C. Delvenne and M. Barahona. Laplacian Dynamics and
 * Multiscale Modular Structure in Networks,
 * doi:10.1109/TNSE.2015.2391998.
 * https://arxiv.org/abs/0812.1770
 *
 * Traag, V. A., et al. « From Louvain to Leiden: Guaranteeing Well-Connected
 * Communities ». Scientific Reports, vol. 9, no 1, décembre 2019, p. 5233.
 * DOI.org (Crossref), doi:10.1038/s41598-019-41695-z.
 * https://arxiv.org/abs/1810.08473
 *
 * [References]
 * https://github.com/vtraag/gephi-leiden-plugin
 * https://github.com/vtraag/leidenalg
 */
var defaults = require('lodash/defaultsDeep'),
    isGraph = require('graphology-utils/is-graph'),
    inferType = require('graphology-utils/infer-type'),
    SparseMap = require('mnemonist/sparse-map'),
    SparseQueueSet = require('mnemonist/sparse-queue-set'),
    createRandomIndex = require('pandemonium/random-index').createRandomIndex,
    utils = require('./utils.js');

var indices = require('graphology-indices/neighborhood/louvain');
var addWeightToCommunity = utils.addWeightToCommunity;

var UndirectedLouvainIndex = indices.UndirectedLouvainIndex,
    DirectedLouvainIndex = indices.DirectedLouvainIndex;

var DEFAULTS = {
  attributes: {
    community: 'community',
    weight: 'weight'
  },
  randomness: 0.01,
  randomWalk: true,
  resolution: 1,
  rng: Math.random,
  weighted: false
};

var EPSILON = 1e-10;

function tieBreaker(bestCommunity, currentCommunity, targetCommunity, delta, bestDelta) {
  if (Math.abs(delta - bestDelta) < EPSILON) {
    if (bestCommunity === currentCommunity) {
      return false;
    }
    else {
      return targetCommunity > bestCommunity;
    }
  }
  else if (delta > bestDelta) {
    return true;
  }

  return false;
}

// TODO: port to graphology-indices?
// TODO: use in zoomOut?
function groupCommunities(index) {
  // To avoid relying on a multimap, we'll use counting sort
  var PointerArray = index.counts.constructor;

  var offsets = new PointerArray(index.C);
  var sorted = new PointerArray(index.C);
  var bounds = new PointerArray(index.C - index.U + 1); // NOTE: could save up

  var n, i, c, b, o;

  n = 0;
  o = 0;

  for (i = 0; i < index.C; i++) {
    c = index.counts[i];

    if (c !== 0) {
      bounds[o++] = n;
      n += c;
      offsets[i] = n;
    }
  }

  bounds[o] = n;

  o = 0;

  for (i = 0; i < index.C; i++) {
    b = index.belongings[i];
    o = --offsets[b];
    sorted[o] = i;
  }

  // var assert = require('assert');
  // var j, l;

  // for (i = 0; i < bounds.length - 1; i++) {
  //   for (j = bounds[i] + 1, l = bounds[i + 1]; j < l; j++) {
  //     assert.strictEqual(index.belongings[sorted[j]], index.belongings[sorted[j - 1]]);
  //   }
  // }

  return [bounds, sorted];
}

function mergeNodesSubset(randomIndex, index, nodes, start, stop) {
  var order = stop - start;
  var resolution = index.resolution;

  var WeightsArray = index.weights.constructor;
  var NodesPointerArray = index.counts.constructor;

  // Counters
  // TODO: consider reusing those from an indexing standpoint
  var totalNodeWeight = 0;
  var clusterWeights = new WeightsArray(order);
  var nonSingletonClusters = new Uint8Array(order);
  var externalEdgeWeightPerCluster = new WeightsArray(order);
  var nodesOrder = new NodesPointerArray(order);
  var belongings = new NodesPointerArray(order);

  var communities = new SparseMap(order);
  var cumulativeIncrement = new Float64Array(order);

  // Iteration variables
  var i, j, n, l, ei, el, et, w;
  var currentCommunity = index.belongings[nodes[start]];

  // Initializing counters
  for (j = 0, i = start; i < stop; i++, j++) {
    n = nodes[i];
    nodesOrder[j] = n;
    belongings[j] = j;
    ei = index.starts[n];
    el = index.starts[n + 1];

    clusterWeights[j] += index.loops[n];

    for(; ei < el; ei++) {
      et = index.neighborhood[ei];

      // Only considering links from the same community
      // TODO: need to recreate a sub neighorhood here...
      if (index.belongings[et] !== currentCommunity)
        continue;

      w = index.weights[et];
      totalNodeWeight += w;
      clusterWeights[j] += w;
      externalEdgeWeightPerCluster[j] += w;
    }
  }

  var s, ri;

  ri = randomIndex(order);

  for (s = 0; s < order; s++, ri++) {
    j = ri % l;
    n = nodesOrder[j];

    // Removing node from its current cluster
    clusterWeights[j] = 0;
    externalEdgeWeightPerCluster[j] = 0;

    // If node is not in a singleton anymore, we can skip
    if (nonSingletonClusters[j] === 1)
      continue;

    // If connectivity constraint is not satisfied, we can skip
    if (externalEdgeWeightPerCluster[j] < clusterWeights[j] * (totalNodeWeight - clusterWeights[j]) * resolution)
      continue

    // Finding neighboring clusters (including the current one)
    communities.clear();
    communities.set(j, 0);

    ei = index.starts[n];
    el = index.starts[n + 1];

    for (; ei < el; ei++) {
      et = index.neighborhood[ei];

      // Only considering links from the same community
      if (index.belongings[et] !== currentCommunity)
        continue;

      // TODO: cannot continue without truncate neighborhood here...
      // TODO: maybe rely on offset from start rather than reindexation
    }
  }
}

function refinePartition(randomIndex, index) {

  // First we need to group by community
  var result = groupCommunities(index);
  var bounds = result[0];
  var nodes = result[1];

  var i, l, start, stop, subpartition;

  for (i = 0, l = bounds.length - 1; i < l; i++) {
    start = bounds[i];
    stop = bounds[i + 1];

    subpartition = mergeNodesSubset(randomIndex, index, nodes, start, stop);
  }
}

function undirectedLeiden(detailed, graph, options) {
  var index = new UndirectedLouvainIndex(graph, {
    attributes: {
      weight: options.attributes.weight
    },
    keepDendrogram: detailed,
    resolution: options.resolution,
    weighted: options.weighted
  });

  var randomIndex = createRandomIndex(options.rng);

  // State variables
  var moveWasMade = true;

  // Communities
  var currentCommunity, targetCommunity;
  var communities = new SparseMap(Float64Array, index.C);

  // Traversal
  var queue = new SparseQueueSet(index.C),
      start,
      end,
      weight,
      ci,
      ri,
      s,
      i,
      j,
      l;

  // Metrics
  var degree,
      targetCommunityDegree;

  // Moves
  var bestCommunity,
      bestDelta,
      deltaIsBetter,
      delta;

  // Details
  var deltaComputations = 0,
      nodesVisited = 0,
      moves = [],
      currentMoves;

  while (moveWasMade) {
    l = index.C;

    moveWasMade = false;
    currentMoves = 0;

    // Traversal of the graph
    ri = options.randomWalk ? randomIndex(l) : 0;

    for (s = 0; s < l; s++, ri++) {
      i = ri % l;
      queue.enqueue(i);
    }

    while (queue.size !== 0) {
      i = queue.dequeue();
      nodesVisited++;

      degree = 0;
      communities.clear();

      currentCommunity = index.belongings[i];

      start = index.starts[i];
      end = index.starts[i + 1];

      // Traversing neighbors
      for (; start < end; start++) {
        j = index.neighborhood[start];
        weight = index.weights[start];

        targetCommunity = index.belongings[j];

        // Incrementing metrics
        degree += weight;
        addWeightToCommunity(communities, targetCommunity, weight);
      }

      // Finding best community to move to
      bestDelta = index.fastDeltaWithOwnCommunity(
        i,
        degree,
        communities.get(currentCommunity) || 0,
        currentCommunity
      );
      bestCommunity = currentCommunity;

      for (ci = 0; ci < communities.size; ci++) {
        targetCommunity = communities.dense[ci];

        if (targetCommunity === currentCommunity)
          continue;

        targetCommunityDegree = communities.vals[ci];

        deltaComputations++;

        delta = index.fastDelta(
          i,
          degree,
          targetCommunityDegree,
          targetCommunity
        );

        deltaIsBetter = tieBreaker(
          bestCommunity,
          currentCommunity,
          targetCommunity,
          delta,
          bestDelta
        );

        if (deltaIsBetter) {
          bestDelta = delta;
          bestCommunity = targetCommunity;
        }
      }

      // Should we move the node?
      if (bestDelta < 0 || bestCommunity !== currentCommunity) {

        // NOTE: this is to allow nodes to move back to their own singleton
        // This code however only deals with modularity (e.g. the condition
        // about bestDelta < 0, which is the delta for moving back to
        // singleton wrt. modularity). Indeed, rarely, the Louvain
        // algorithm can produce such cases when a node would be better in
        // a singleton that in its own community when considering self loops
        // or a resolution != 1. In this case, delta with your own community
        // is indeed less than 0. To handle different metrics, one should
        // consider computing the delta for going back to singleton because
        // it might not be 0.
        if (bestDelta < 0) {
          index.isolate(i, degree);
        }
        else {
          index.move(i, degree, bestCommunity);
        }

        moveWasMade = true;
        currentMoves++;

        // Adding neighbors from other communities to the queue
        start = index.starts[i];
        end = index.starts[i + 1];

        for (; start < end; start++) {
          j = index.neighborhood[start];
          targetCommunity = index.belongings[j];

          if (targetCommunity !== bestCommunity)
            queue.enqueue(j);
        }
      }
    }

    moves.push(currentMoves);

    // We continue working on the induced graph
    if (moveWasMade) {
      console.time('refine');
      refinePartition(randomIndex, index);
      console.timeEnd('refine');
      throw new Error('unimplemented');
      index.zoomOut();
    }
  }

  var results = {
    index: index,
    deltaComputations: deltaComputations,
    nodesVisited: nodesVisited,
    moves: moves
  };

  return results;
}

function directedLeiden(detailed, graph, options) {
  var index = new DirectedLouvainIndex(graph, {
    attributes: {
      weight: options.attributes.weight
    },
    keepDendrogram: detailed,
    resolution: options.resolution,
    weighted: options.weighted
  });

  var randomIndex = createRandomIndex(options.rng);

  // State variables
  var moveWasMade = true;

  // Communities
  var currentCommunity, targetCommunity;
  var communities = new SparseMap(Float64Array, index.C);

  // Traversal
  var queue = new SparseQueueSet(index.C),
      start,
      end,
      offset,
      out,
      weight,
      ci,
      ri,
      s,
      i,
      j,
      l;

  // Metrics
  var inDegree,
      outDegree,
      targetCommunityDegree;

  // Moves
  var bestCommunity,
      bestDelta,
      deltaIsBetter,
      delta;

  // Details
  var deltaComputations = 0,
      nodesVisited = 0,
      moves = [],
      currentMoves;

  while (moveWasMade) {
    l = index.C;

    moveWasMade = false;
    currentMoves = 0;

    // Traversal of the graph
    ri = options.randomWalk ? randomIndex(l) : 0;

    for (s = 0; s < l; s++, ri++) {
      i = ri % l;
      queue.enqueue(i);
    }

    while (queue.size !== 0) {
      i = queue.dequeue();
      nodesVisited++;

      inDegree = 0;
      outDegree = 0;
      communities.clear();

      currentCommunity = index.belongings[i];

      start = index.starts[i];
      end = index.starts[i + 1];
      offset = index.offsets[i];

      // Traversing neighbors
      for (; start < end; start++) {
        out = start < offset;
        j = index.neighborhood[start];
        weight = index.weights[start];

        targetCommunity = index.belongings[j];

        // Incrementing metrics
        if (out)
          outDegree += weight;
        else
          inDegree += weight;

        addWeightToCommunity(communities, targetCommunity, weight);
      }

      // Finding best community to move to
      bestDelta = index.deltaWithOwnCommunity(
        i,
        inDegree,
        outDegree,
        communities.get(currentCommunity) || 0,
        currentCommunity
      );
      bestCommunity = currentCommunity;

      for (ci = 0; ci < communities.size; ci++) {
        targetCommunity = communities.dense[ci];

        if (targetCommunity === currentCommunity)
          continue;

        targetCommunityDegree = communities.vals[ci];

        deltaComputations++;

        delta = index.delta(
          i,
          inDegree,
          outDegree,
          targetCommunityDegree,
          targetCommunity
        );

        deltaIsBetter = tieBreaker(
          bestCommunity,
          currentCommunity,
          targetCommunity,
          delta,
          bestDelta
        );

        if (deltaIsBetter) {
          bestDelta = delta;
          bestCommunity = targetCommunity;
        }
      }

      // Should we move the node?
      if (bestDelta < 0 || bestCommunity !== currentCommunity) {

        // NOTE: this is to allow nodes to move back to their own singleton
        // This code however only deals with modularity (e.g. the condition
        // about bestDelta < 0, which is the delta for moving back to
        // singleton wrt. modularity). Indeed, rarely, the Louvain
        // algorithm can produce such cases when a node would be better in
        // a singleton that in its own community when considering self loops
        // or a resolution != 1. In this case, delta with your own community
        // is indeed less than 0. To handle different metrics, one should
        // consider computing the delta for going back to singleton because
        // it might not be 0.
        if (bestDelta < 0) {
          index.isolate(i, inDegree, outDegree);
        }
        else {
          index.move(i, inDegree, outDegree, bestCommunity);
        }

        moveWasMade = true;
        currentMoves++;

        // Adding neighbors from other communities to the queue
        start = index.starts[i];
        end = index.starts[i + 1];

        for (; start < end; start++) {
          j = index.neighborhood[start];
          targetCommunity = index.belongings[j];

          if (targetCommunity !== bestCommunity)
            queue.enqueue(j);
        }
      }
    }

    moves.push(currentMoves);

    // We continue working on the induced graph
    if (moveWasMade)
      index.zoomOut();
  }

  var results = {
    index: index,
    deltaComputations: deltaComputations,
    nodesVisited: nodesVisited,
    moves: moves
  };

  return results;
}

/**
 * Function returning the communities mapping of the graph.
 *
 * @param  {boolean} assign             - Assign communities to nodes attributes?
 * @param  {boolean} detailed           - Whether to return detailed information.
 * @param  {Graph}   graph              - Target graph.
 * @param  {object}  options            - Options:
 * @param  {object}    attributes         - Attribute names:
 * @param  {string}      community          - Community node attribute name.
 * @param  {string}      weight             - Weight edge attribute name.
 * @param  {string}    deltaComputation   - Method to use to compute delta computations.
 * @param  {number}    randomness         - Randomness parameter.
 * @param  {boolean}   randomWalk         - Whether to traverse the graph in random order.
 * @param  {number}    resolution         - Resolution parameter.
 * @param  {function}  rng                - RNG function to use.
 * @param  {boolean}   weighted           - Whether to compute the weighted version.
 * @return {object}
 */
function leiden(assign, detailed, graph, options) {
  if (!isGraph(graph))
    throw new Error('graphology-communities-leiden: the given graph is not a valid graphology instance.');

  var type = inferType(graph);

  if (type === 'mixed')
    throw new Error('graphology-communities-leiden: cannot run the algorithm on a true mixed graph.');

  // Attributes name
  options = defaults({}, options, DEFAULTS);

  // Empty graph case
  var c = 0;

  if (graph.size === 0) {
    if (assign) {
      graph.forEachNode(function(node) {
        graph.setNodeAttribute(node, options.attributes.communities, c++);
      });

      return;
    }

    var communities = {};

    graph.forEachNode(function(node) {
      communities[node] = c++;
    });

    if (!detailed)
      return communities;

    return {
      communities: communities,
      count: graph.order,
      deltaComputations: 0,
      dendrogram: null,
      level: 0,
      modularity: NaN,
      moves: null,
      nodesVisited: 0,
      resolution: options.resolution
    };
  }

  var fn = type === 'undirected' ? undirectedLeiden : directedLeiden;

  var results = fn(detailed, graph, options);

  var index = results.index;

  // Standard output
  if (!detailed) {
    if (assign) {
      index.assign(options.attributes.community);
      return;
    }

    return index.collect();
  }

  // Detailed output
  var output = {
    count: index.C,
    deltaComputations: results.deltaComputations,
    dendrogram: index.dendrogram,
    level: index.level,
    modularity: index.modularity(),
    moves: results.moves,
    nodesVisited: results.nodesVisited,
    resolution: options.resolution
  };

  if (assign) {
    index.assign(options.attributes.community);
    return output;
  }

  output.communities = index.collect();

  return output;
}

/**
 * Exporting.
 */
var fn = leiden.bind(null, false, false);
fn.assign = leiden.bind(null, true, false);
fn.detailed = leiden.bind(null, false, true);
fn.defaults = DEFAULTS;

module.exports = fn;
