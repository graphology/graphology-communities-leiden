/**
 * Graphology Leiden Utils
 * ========================
 *
 * Miscellaneous utilities used by the Leiden algorithm.
 */
function UndirectedLeidenAddenda(index, options) {
  options = options || {};

  this.index = index;

  var NodesPointerArray = index.counts.constructor;

  // Used to group nodes by communities
  this.communitiesOffsets = new NodesPointerArray(index.C);
  this.nodesSortedByCommunities = new NodesPointerArray(index.C);
  this.communitiesBounds = new NodesPointerArray(index.C + 1);
  this.B = 0;
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

exports.UndirectedLeidenAddenda = UndirectedLeidenAddenda;
