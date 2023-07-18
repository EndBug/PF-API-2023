#include <limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BLACK false
#define RED true

#define FORWARD false
#define BACKWARD true

#define INFINITY INT_MAX - 10

#define DEBUG 0

/** Heap structure that contains the vehicle info */
typedef struct vehicles_ {
  /** Number of elements in the heap */
  int size;

  /** Array representation of the heap */
  int *heap;
} vehicles_t;

/** Node of the RB tree that represents a station */
typedef struct station_ {
  /** Distance of the station form the start of the road */
  int key;

  /** Color of the station in the RB tree */
  bool color;

  /** Pointer to the left child of the station in the RB tree */
  struct station_ *left;

  /** Pointer to the right child of the station in the RB tree */
  struct station_ *right;

  /** Pointer to the parent of the station in the RB tree */
  struct station_ *parent;

  /** Pointer to the vehicles heap structure */
  vehicles_t *vehicles;
} station_t;

/** RB tree that contains all stations */
typedef struct station_tree_ {
  /** Root of the RB tree */
  station_t *root;

  /** Special node used as the parent of the root and the child of all leaves */
  station_t *nil;

  /** Number of nodes in the RB tree */
  int size;
} station_tree_t;

/** Node of the graph used by the pathfinding algorithm*/
typedef struct graph_node_ {
  /** Distance of the station from the start of the road */
  int key;

  /** Maximum range of the vehicles in the station */
  int range;

  /** Distance from the source node */
  int d;

  /** Index of the previous node in the path */
  int p;
} graph_node_t;

/** Directed Acyclic Graph that represents the road */
typedef struct graph_ {
  /** Array of pointers to each graph node */
  graph_node_t **nodes;

  /** Number of nodes in the graph */
  int size;
} graph_t;

void cmd_add_station(station_tree_t *T);
void cmd_remove_station(station_tree_t *T);
void cmd_add_vehicle(station_tree_t *T);
void cmd_remove_vehicle(station_tree_t *T);
void cmd_plan_trip(station_tree_t *T);
void print_stations(station_tree_t *T, int v);
void print_rb_tree(station_tree_t *T);
void destroy_tree(station_tree_t *T);

int main() {
  station_tree_t *T = malloc(sizeof(station_tree_t));

  if (T) {
    T->nil = malloc(sizeof(station_t));
    T->nil->left = T->nil;
    T->nil->right = T->nil;
    T->nil->parent = T->nil;
    T->nil->color = BLACK;
    T->root = T->nil;
    T->size = 0;

    char command[20];
    // eof_res should be -1 when there is no stdin left to read
    int eof_res = scanf("%s", command);

    while (eof_res >= 0) {
      if (!strcmp(command, "aggiungi-stazione"))
        cmd_add_station(T);
      else if (!strcmp(command, "demolisci-stazione"))
        cmd_remove_station(T);
      else if (!strcmp(command, "aggiungi-auto"))
        cmd_add_vehicle(T);
      else if (!strcmp(command, "rottama-auto"))
        cmd_remove_vehicle(T);
      else if (!strcmp(command, "pianifica-percorso"))
        cmd_plan_trip(T);
      else {
        printf("Unknown command: %s\n", command);
        return -1;
      }

      eof_res = scanf("%s", command);
    }

    destroy_tree(T);
  } else {
    printf("[main] Allocation error.\n");
    return -1;
  }
}

// #region max_heap
/** Returns the 1-based index of the left child of the n-th element of a 1-based
 * array */
int left(int n) { return 2 * n; }
/** Returns the 1-based index of the right child of the n-th element of a
 * 1-based array */
int right(int n) { return 2 * n + 1; }
/** Returns the 1-based index of the parent of the n-th element of a 1-based
 * array */
int parent(int n) { return n / 2; }

/** Swaps two integer values */
void swap(int *a, int *b) {
  int tmp = *a;
  *a = *b;
  *b = tmp;
}

/**
 * @brief Makes a value descend the heap towards the leaves until it's greater
 * than the children
 *
 * @param A A pointer to the heap structure
 * @param n The 1-based index of the element in the heap
 */
void max_heapify(vehicles_t *A, int n) {
  int l = left(n);
  int r = right(n);
  int posmax;

  if (l <= A->size && A->heap[l - 1] > A->heap[n - 1])
    posmax = l;
  else
    posmax = n;

  if (r <= A->size && A->heap[r - 1] > A->heap[posmax - 1])
    posmax = r;

  if (posmax != n) {
    swap(&A->heap[n - 1], &A->heap[posmax - 1]);
    max_heapify(A, posmax);
  }
}

/**
 * @brief Transforms an array into the representation of a max heap
 *
 * @param A A pointer to the heap structure
 */
void build_max_heap(vehicles_t *A) {
  int i;
  for (i = A->size / 2; i > 0; i--) // i is 1-based!
    max_heapify(A, i);
}

/**
 * @brief Removes the maximum element of a heap and returns it
 *
 * @param A A pointer to the heap structure
 * @return The maximum element of the heap, or -1 if the heap is empty
 */
int heap_shift(vehicles_t *A) {
  if (A->size < 1)
    return -1;

  int max = A->heap[0];
  A->heap[0] = A->heap[A->size - 1];
  A->size = A->size - 1;

  max_heapify(A, 1);

  return max;
}

/**
 * @brief Add an element to the heap
 *
 * @param A A pointer to the heap structure
 * @param key The key of the new element
 */
void heap_insert(vehicles_t *A, int key) {
  A->size = A->size + 1;
  A->heap = realloc(A->heap, sizeof(int) * A->size);
  A->heap[A->size - 1] = key;

  int i = A->size;
  while (i > 1 && A->heap[parent(i) - 1] < A->heap[i - 1]) {
    swap(&A->heap[parent(i) - 1], &A->heap[i - 1]);
    i = parent(i);
  }
}

/**
 * @brief Returns the maximum element of the heap
 *
 * @param A A pointer to the heap structure
 * @return The maximum element, or -1 if the heap is empty
 */
int heap_max(vehicles_t *A) {
  if (A->size < 1)
    return -1;

  return A->heap[0];
}

/**
 * @brief Deletes an element from the heap
 *
 * @param A A pointer to the heap structure
 * @param key The value f the element to delete
 * @return Whether the element was deleted
 */
bool heap_delete(vehicles_t *A, int key) {
  int i;
  for (i = 0; i < A->size && A->heap[i] != key; i++)
    ;

  if (i == A->size)
    return false;

  A->heap[i] = A->heap[A->size - 1];
  A->size = A->size - 1;

  max_heapify(A, i + 1); // 1-based

  return true;
}

// #endregion

// #region RB_tree

/**
 * @brief Creates a station
 *
 * @param T A pointer to the tree the station will be inserted in (please note
 * that this function will not do it, you need to use `insert_station`)
 * @param key The key of the station (its distance from the start of the road)
 * @param vehicles An array of vehicles to add
 * @param vn The number of vehicles in the array
 * @return A pointer to the station node
 */
station_t *create_station(station_tree_t *T, int key, int *vehicles, int vn) {
  station_t *s = malloc(sizeof(station_t));

  if (s) {
    s->key = key;
    s->color = RED;
    s->left = T->nil;
    s->right = T->nil;
    s->parent = T->nil;

    s->vehicles = malloc(sizeof(vehicles_t));
    if (s->vehicles) {
      s->vehicles->heap = vehicles;
      s->vehicles->size = vn;
      build_max_heap(s->vehicles);
    } else
      printf("[create_station] Allocation error.\n");
  } else
    printf("[create_station] Allocation error.\n");

  return s;
}

/**
 * @brief Rotates an RB tree to the left
 *
 * @param T A pointer to the station tree to rotate
 * @param x A pointer to the pivot of the rotation
 */
void station_rotate_left(station_tree_t *T, station_t *x) {
  station_t *y = x->right;
  x->right = y->left;

  if (y->left != T->nil)
    y->left->parent = x;
  y->parent = x->parent;

  if (x->parent == T->nil)
    T->root = y;
  else if (x == x->parent->left)
    x->parent->left = y;
  else
    x->parent->right = y;

  y->left = x;
  x->parent = y;
}

/**
 * @brief Rotates an RB tree to the right
 *
 * @param T A pointer to the station tree to rotate
 * @param x A pointer to the pivot of the rotation
 */
void station_rotate_right(station_tree_t *T, station_t *x) {
  station_t *y = x->left;
  x->left = y->right;

  if (y->right != T->nil)
    y->right->parent = x;
  y->parent = x->parent;

  if (x->parent == T->nil)
    T->root = y;
  else if (x == x->parent->right)
    x->parent->right = y;
  else
    x->parent->left = y;

  y->right = x;
  x->parent = y;
}

/**
 * @brief Fixes the station RB tree after insertion
 *
 * @param T A pointer to the stations tree
 * @param z A pointer to the node that has just been inserted via
 * `insert_station`
 */
void insert_station_fixup(station_tree_t *T, station_t *z) {
  station_t *y;

  while (z->parent->color == RED) {
    if (z->parent == z->parent->parent->left) {
      y = z->parent->parent->right;

      if (y->color == RED) {
        z->parent->color = BLACK;
        y->color = BLACK;
        z->parent->parent->color = RED;
        z = z->parent->parent;
      } else {
        if (z == z->parent->right) {
          z = z->parent;
          station_rotate_left(T, z);
        }

        z->parent->color = BLACK;
        z->parent->parent->color = RED;
        station_rotate_right(T, z->parent->parent);
      }
    } else {
      y = z->parent->parent->left;

      if (y->color == RED) {
        z->parent->color = BLACK;
        y->color = BLACK;
        z->parent->parent->color = RED;
        z = z->parent->parent;
      } else {
        if (z == z->parent->left) {
          z = z->parent;
          station_rotate_right(T, z);
        }

        z->parent->color = BLACK;
        z->parent->parent->color = RED;
        station_rotate_left(T, z->parent->parent);
      }
    }
  }

  T->root->color = BLACK;
}

/**
 * Inserts a station into the RB tree
 * @param *T A pointer to a station RB tree
 * @param *x A pointer to the station node to insert, generated by
 * `create_station`
 */
void insert_station(station_tree_t *T, station_t *z) {
  station_t *y = T->nil;
  station_t *x = T->root;

  while (x != T->nil) {
    y = x;
    if (z->key < x->key)
      x = x->left;
    else
      x = x->right;
  }

  z->parent = y;

  if (y == T->nil)
    T->root = z;
  else if (z->key < y->key)
    y->left = z;
  else
    y->right = z;

  z->left = T->nil;
  z->right = T->nil;
  z->color = RED;

  T->size = T->size + 1;
  insert_station_fixup(T, z);
}

/**
 * @brief Finds a station node in the tree
 *
 * @param T A pointer to the tree to search
 * @param key The key of the node to search for
 * @return A pointer to the target node, or to T->nil when it can't be found
 */
station_t *find_station(station_tree_t *T, int key) {
  station_t *curr = T->root;

  while (curr != T->nil && curr->key != key) {
    if (curr->key < key)
      curr = curr->right;
    else
      curr = curr->left;
  }

  return curr;
}

/** Recursive step of `print_stations` */
void print_stations_rec(station_tree_t *T, station_t *curr, int v) {
  if (curr->left != T->nil)
    print_stations_rec(T, curr->left, v);

  printf("%d: ", curr->key);
  if (v) {
    printf("\n  ");
    int i;
    for (i = 0; curr->vehicles != NULL && i < curr->vehicles->size; i++)
      printf("%d ", curr->vehicles->heap[i]);
    printf("\n");
  }

  if (curr->right != T->nil)
    print_stations_rec(T, curr->right, v);
}
/** Prints all the stations in the RB tree */
void print_stations(station_tree_t *T, int full) {
  return print_stations_rec(T, T->root, full);
}

void print_rb_tree_rec(station_tree_t *T, station_t *node, int level) {
  char col;
  int i;
  level++;

  if (node->color == BLACK)
    col = 'B';
  else
    col = 'R';
  printf("--%d%c", node->key, col);

  if (node->right != T->nil)
    print_rb_tree_rec(T, node->right, level);
  if (node->left != T->nil) {
    printf("\n");
    for (i = 0; i < level; i++)
      printf("      ");
    print_rb_tree_rec(T, node->left, level);
  }
  // printf("\n");
}
void print_rb_tree(station_tree_t *T) { print_rb_tree_rec(T, T->root, 0); }

/** Subroutine called by `remove_station` */
void station_transplant(station_tree_t *T, station_t *u, station_t *v) {
  if (u->parent == T->nil)
    T->root = v;
  else if (u == u->parent->left)
    u->parent->left = v;
  else
    u->parent->right = v;

  v->parent = u->parent;
}

/**
 * @brief Finds the element of a subtree with the lowest key
 *
 * @param T A pointer to the stations tree
 * @param x A pointer to the root of the subtree to consider
 * @return A pointer to the element with the lowest key, the leftmost child
 */
station_t *tree_minimum(station_tree_t *T, station_t *x) {
  while (x->left != T->nil)
    x = x->left;
  return x;
}

/**
 * @brief Fixes the RB station tree after deletion
 *
 * @param T A pointer to the stations tree
 * @param x A pointer to the removed node
 */
void remove_station_fixup(station_tree_t *T, station_t *x) {
  station_t *w;

  while (x != T->root && x->color == BLACK) {
    if (x == x->parent->left) {
      w = x->parent->right;

      if (w->color == RED) {
        w->color = BLACK;
        x->parent->color = RED;
        station_rotate_left(T, x->parent);
        w = x->parent->right;
      }

      if (w->left->color == BLACK && w->right->color == BLACK) {
        w->color = RED;
        x = x->parent;
      } else {
        if (w->right->color == BLACK) {
          w->left->color = BLACK;
          w->color = RED;
          station_rotate_right(T, w);
          w = x->parent->right;
        }

        w->color = x->parent->color;
        x->parent->color = BLACK;
        w->right->color = BLACK;
        station_rotate_left(T, x->parent);
        x = T->root;
      }
    } else {
      w = x->parent->left;

      if (w->color == RED) {
        w->color = BLACK;
        x->parent->color = RED;
        station_rotate_right(T, x->parent);
        w = x->parent->left;
      }

      if (w->right->color == BLACK && w->left->color == BLACK) {
        w->color = RED;
        x = x->parent;
      } else {
        if (w->left->color == BLACK) {
          w->right->color = BLACK;
          w->color = RED;
          station_rotate_left(T, w);
          w = x->parent->left;
        }

        w->color = x->parent->color;
        x->parent->color = BLACK;
        w->left->color = BLACK;
        station_rotate_right(T, x->parent);
        x = T->root;
      }
    }
  }

  x->color = BLACK;
}

/**
 * @brief Deallocates a station
 *
 * @param z A pointer to the station to deallocate
 */
void destroy_station(station_t *z) {
  free(z->vehicles->heap);
  free(z->vehicles);
  free(z);
}

/**
 * @brief Removes a station form the RB tree
 *
 * @param T A pointer to the stations tree
 * @param z A pointer to the node to remove
 */
void remove_station(station_tree_t *T, station_t *z) {
  station_t *x;
  station_t *y = z;
  int y_orig_col = y->color;

  if (z->left == T->nil) {
    x = z->right;
    station_transplant(T, z, z->right);
  } else if (z->right == T->nil) {
    x = z->left;
    station_transplant(T, z, z->left);
  } else {
    y = tree_minimum(T, z->right);
    y_orig_col = y->color;
    x = y->right;

    if (y->parent == z)
      x->parent = y;
    else {
      station_transplant(T, y, y->right);
      y->right = z->right;
      y->right->parent = y;
    }

    station_transplant(T, z, y);
    y->left = z->left;
    y->left->parent = y;
    y->color = z->color;
  }

  T->size = T->size - 1;
  if (y_orig_col == BLACK)
    remove_station_fixup(T, x);
}

/**
 * @brief Deallocates stations recursively
 *
 * @param T A pointer to the staations tree
 * @param curr The current station node
 */
void destroy_tree_rec(station_tree_t *T, station_t *curr) {
  if (curr != T->nil) {
    destroy_tree_rec(T, curr->left);
    destroy_tree_rec(T, curr->right);

    destroy_station(curr);
  }
}
/**
 * @brief Deallocates a stations tree from memory
 *
 * @param T A pointer to the station tree
 */
void destroy_tree(station_tree_t *T) {
  destroy_tree_rec(T, T->root);
  free(T->nil);
  free(T);
}

// #endregion

// #region graph
/**
 * @brief Get the direction of a graph
 *
 * @param source The distance of the source station from the start of the road
 * @param dest The distance of the destination station from the start of the
 * road
 * @return Either `FORWARD` or `BACKWARD`
 */
bool get_direction(int source, int dest) {
  if (source <= dest)
    return FORWARD;
  else
    return BACKWARD;
}

/**
 * @brief Initializes the source of the trip for `dag_shortest_paths`
 *
 * @param G A pointer to the graph
 * @param s A pointer to the source node
 */
void initialize_single_source(graph_t *G) {
  int i;
  for (i = 0; i < G->size; i++) {
    G->nodes[i]->d = INFINITY;
    G->nodes[i]->p = -1;
  }

  G->nodes[0]->d = 0;
}

/**
 * @brief Performs relaxation on two given graph nodes for `dag_shortest_paths`
 *
 * @param u A pointer to the local source node
 * @param v A pointer to the local destination node
 * @param direction The direction in which the road is being used
 */
void relax(graph_t *G, int ui, int vi, bool direction) {
  // We're considering all weights to be 1

  graph_node_t *u = G->nodes[ui];
  graph_node_t *v = G->nodes[vi];

  // When ==, choose the smaller key: this ensures that the algorithm always
  // chooses the stations closest to the start of the road, even when going
  // backwards.
  if (v->d > u->d + 1 || (v->d == u->d + 1 && G->nodes[v->p]->key > u->key)) {
    if (DEBUG)
      printf("relax %d(%d) %d(%d)\n", u->key, u->d, v->key, v->d);

    v->d = u->d + 1;
    v->p = ui;
  }
}

/**
 * Algorithm to compute the shortest path in a Directed Acyclic Graph
 * The `relax` procedure has been modified to prefer nodes closer to the start
 * of the road, as per specification.
 *
 * @param G A pointer to the graph
 * @param direction The direction in which the road is being used
 */
void dag_shortest_path(graph_t *G, bool direction) {
  initialize_single_source(G);

  int ui;
  for (ui = 0; ui < G->size; ui++) {
    graph_node_t *u = G->nodes[ui];

    int vi;
    if (direction == FORWARD) {
      for (vi = ui + 1; vi < G->size && G->nodes[vi]->key <= u->key + u->range;
           vi++)
        relax(G, ui, vi, direction);
    } else {
      for (vi = ui + 1; vi < G->size && G->nodes[vi]->key >= u->key - u->range;
           vi++)
        relax(G, ui, vi, direction);
    }
  }
}

/**
 * Recursively adds nodes to the graph
 *
 * @param G A pointer to the graph
 * @param T A pointer to the stations tree
 * @param curr A pointer to the current tree node
 * @param source The distance of the source station from the start of the road
 * @param dest The distance of the destination station from the start of the
 * road
 * @param direction The direction in which the road is being used
 */
void add_nodes(graph_t *G, station_tree_t *T, station_t *curr, int source,
               int dest, bool direction) {
  if (curr == T->nil)
    return;

  if (direction == FORWARD) {
    if (curr->key > source)
      add_nodes(G, T, curr->left, source, dest, direction);

    if (curr->key >= source && curr->key <= dest) {
      G->nodes[G->size] = malloc(sizeof(graph_node_t));

      if (G->nodes[G->size]) {
        G->nodes[G->size]->key = curr->key;
        G->nodes[G->size]->range = heap_max(curr->vehicles);
        G->size = G->size + 1;
      } else
        printf("[add_nodes] Allocation error.\n");
    }

    if (curr->key < dest)
      add_nodes(G, T, curr->right, source, dest, direction);
  } else {
    if (curr->key < source)
      add_nodes(G, T, curr->right, source, dest, direction);

    if (curr->key >= dest && curr->key <= source) {
      G->nodes[G->size] = malloc(sizeof(graph_node_t));

      if (G->nodes[G->size]) {
        G->nodes[G->size]->key = curr->key;
        G->nodes[G->size]->range = heap_max(curr->vehicles);
        G->size = G->size + 1;
      } else
        printf("[add_nodes] Allocation error.\n");
    }

    if (curr->key > dest)
      add_nodes(G, T, curr->left, source, dest, direction);
  }
}

/**
 * Builds a Directed Acyclic Graph from the road
 *
 * @param T A pointer to the stations tree
 * @param source The distance of the source station from the start of the road
 * @param dest The distance of the destination station from the start of the
 * road
 * @param direction The direction in which the road is being used
 * @return A pointer to the new graph
 */
graph_t *build_graph(station_tree_t *T, int source, int dest, bool direction) {
  graph_t *G = malloc(sizeof(graph_t));

  if (G != NULL) {
    G->nodes = malloc(sizeof(graph_node_t) * T->size);
    G->size = 0;

    if (G->nodes != NULL) {
      add_nodes(G, T, T->root, source, dest, direction);

      G->nodes = realloc(G->nodes, sizeof(graph_node_t) * G->size);
    } else
      printf("[build_graph] Allocation error.\n");
  } else
    printf("[build_graph] Allocation error.\n");

  return G;
}

/** Deallocates a graph from memory */
void destroy_graph(graph_t *G) {
  int i;
  for (i = 0; i < G->size; i++)
    free(G->nodes[i]);

  free(G->nodes);
  free(G);
}
// #endregion

// #region commands

void cmd_add_station(station_tree_t *T) {
  int dist, vn;
  if (scanf("%d %d", &dist, &vn))
    ;

  if (find_station(T, dist) != T->nil) {
    printf("non aggiunta\n");
    char _[50000] = "";
    if (scanf("%[^\n]", _))
      ; // Discards the rest of the line
  } else {
    int i;
    int *vehicles = malloc(sizeof(int) * vn);

    if (vehicles) {
      for (i = 0; i < vn; i++)
        if (scanf("%d", &vehicles[i]))
          ;

      station_t *s = create_station(T, dist, vehicles, vn);

      insert_station(T, s);

      printf("aggiunta\n");
    } else
      printf("[cmd_add_station] Allocation error.\n");
  }
}

void cmd_remove_station(station_tree_t *T) {
  int dist;
  if (scanf("%d", &dist))
    ;

  station_t *node = find_station(T, dist);
  if (node != T->nil) {
    remove_station(T, node);

    destroy_station(node);

    printf("demolita\n");
  } else {
    printf("non demolita\n");
  }
}

void cmd_add_vehicle(station_tree_t *T) {
  int dist, range;
  if (scanf("%d %d", &dist, &range))
    ;

  station_t *station = find_station(T, dist);
  if (station != T->nil) {
    heap_insert(station->vehicles, range);
    printf("aggiunta\n");
  } else
    printf("non aggiunta\n");
}

void cmd_remove_vehicle(station_tree_t *T) {
  int dist, range;
  if (scanf("%d %d", &dist, &range))
    ;

  station_t *station = find_station(T, dist);
  if (station != T->nil && heap_delete(station->vehicles, range))
    printf("rottamata\n");
  else
    printf("non rottamata\n");
}

void cmd_plan_trip(station_tree_t *T) {
  int source, dest;
  if (scanf("%d %d", &source, &dest))
    ;

  bool direction = get_direction(source, dest);
  graph_t *G = build_graph(T, source, dest, direction);
  if (G != NULL && G->size != 0) {
    dag_shortest_path(G, direction);

    if (DEBUG)
      print_stations(T, 1);

    if (G->nodes[G->size - 1]->d != INFINITY) {
      int inv_path[G->size], path_size = 0;
      int last = -1;

      int curr = G->size - 1;
      while (curr != -1) {
        inv_path[path_size] = G->nodes[curr]->key;
        path_size++;

        last = curr;
        curr = G->nodes[curr]->p;
      }

      if (last == -1 || G->nodes[last]->key != source)
        printf("nessun percorso\n");
      else {
        for (; path_size > 0; path_size--) {
          printf("%d", inv_path[path_size - 1]);
          if (path_size > 1)
            printf(" ");
        }
        printf("\n");
      }
    } else
      printf("nessun percorso\n");

    destroy_graph(G);
  } else
    printf("nessun percorso\n");
}

// #endregion
