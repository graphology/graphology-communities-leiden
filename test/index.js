/**
 * Graphology Leiden Unit Tests
 * =============================
 */
var assert = require('chai').assert;
var Graph = require('graphology');
var louvainIndices = require('graphology-indices/neighborhood/louvain');
var leidenIndices = require('../utils.js');
var mergeClique = require('graphology-utils/merge-clique');

var UndirectedLouvainIndex = louvainIndices.UndirectedLouvainIndex;
var UndirectedLeidenAddenda = leidenIndices.UndirectedLeidenAddenda;

// var ARCTIC = require('./resources/arctic.json');

function getDoubleCliqueGraph() {
  var graph = new Graph.UndirectedGraph();
  mergeClique(graph, [0, 1, 2]);
  mergeClique(graph, [3, 4, 5]);
  graph.addEdge(2, 4);

  return graph;
}

describe('graphology-communities-leiden', function() {
  describe('UndirectedLeidenAddenda', function() {
    it('should properly group by communities.', function() {
      var graph = getDoubleCliqueGraph();

      var index = new UndirectedLouvainIndex(graph);
      var addenda = new UndirectedLeidenAddenda(index);

      index.expensiveMove(1, 0);
      index.expensiveMove(2, 0);
      index.expensiveMove(3, 4);
      index.expensiveMove(5, 4);

      addenda.groupByCommunities();

      assert.strictEqual(addenda.B, 2);

      assert.deepStrictEqual(addenda.communities(), [[0, 1, 2], [3, 4, 5]]);
    });

    it('should properly refine the partition.', function() {
      var graph = getDoubleCliqueGraph();

      var index = new UndirectedLouvainIndex(graph);
      var addenda = new UndirectedLeidenAddenda(index);

      index.expensiveMove(1, 0);
      index.expensiveMove(2, 0);
      index.expensiveMove(3, 4);
      index.expensiveMove(5, 4);

      addenda.refinePartition();

      assert.deepStrictEqual(addenda.belongings, new Uint8Array([0, 1, 2, 3, 4, 5]));
    });

    it.skip('should work with fig C1.', function() {

    });
  });
});
