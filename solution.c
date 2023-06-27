#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BLACK 0
#define RED 1

typedef struct vehicles_ {
  int size;
  int *heap;
} vehicles_t;

typedef struct station_ {
  int key, color;
  struct station_ *left, *right, *parent;

  vehicles_t *vehicles;
} station_t;
typedef struct station_tree_ {
  station_t *root, *nil;
} station_tree_t;

void cmd_add_station(station_tree_t *T);
void print_stations(station_tree_t *T);

int main() {
  station_tree_t *T = malloc(sizeof(station_tree_t));

  if (T) {
    T->nil = malloc(sizeof(station_t));
    T->nil->left = T->nil;
    T->nil->right = T->nil;
    T->nil->parent = T->nil;
    T->root = T->nil;

    char command[20];
    // eof_res should be -1 when there is no stdin left to read
    int eof_res = scanf("%s", command);

    while (eof_res >= 0) {
      if (!strcmp(command, "aggiungi-stazione"))
        cmd_add_station(T);
      else if (!strcmp(command, "demolisci-stazione")) {
      } else if (!strcmp(command, "aggiungi-auto")) {
      } else if (!strcmp(command, "rottama-auto")) {
      } else if (!strcmp(command, "pianifica-percorso")) {
      } else {
        printf("Unknown command: %s\n", command);
        return -1;
      }

      eof_res = scanf("%s", command);
    }
  } else {
    printf("[main] Allocation error.\n");
    return -1;
  }
}

#pragma region max_heap
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
 * @return The maximum element of the heap
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
  A->heap[A->size - 1] = key;

  int i = A->size;
  while (i > 1 && A->heap[parent(i) - 1] < A->heap[i - 1]) {
    swap(&A->heap[parent(i) - 1], &A->heap[i - 1]);
    i = parent(i);
  }
}

#pragma endregion

#pragma region RB_tree

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

  while (z->parent != T->nil && z->parent->color == RED) {
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
void insert_station(station_tree_t *T, station_t *x) {
  station_t *pre = T->nil;
  station_t *cur = T->root;

  while (cur != T->nil) {
    pre = cur;
    if (x->key < cur->key)
      cur = cur->left;
    else
      cur = cur->right;
  }

  x->parent = pre;
  if (pre == T->nil)
    T->root = x;
  else if (x->key < pre->key)
    pre->left = x;
  else
    pre->right = x;

  x->color = RED;
  insert_station_fixup(T, x);
}

/**
 * @brief Recursive step of `find_station`
 *
 * @param T A pointer to the tree to search
 * @param key The key of the node to search for
 * @param curr A pointer to the current node
 * @return A pointer to the target node
 */
station_t *find_station_rec(station_tree_t *T, int key, station_t *curr) {
  if (curr == T->nil || curr->key == key)
    return curr;

  station_t *res = T->nil;
  if (curr->left != T->nil)
    res = find_station_rec(T, key, curr->left);
  if (res == T->nil && curr->right != T->nil)
    res = find_station_rec(T, key, curr->right);

  return res;
}
/**
 * @brief Finds a station node in the tree
 *
 * @param T A pointer to the tree to search
 * @param key The key of the node to search for
 * @return A pointer to the target node
 */
station_t *find_station(station_tree_t *T, int key) {
  return find_station_rec(T, key, T->root);
}

void print_stations_rec(station_tree_t *T, station_t *curr) {
  if (curr->left != T->nil)
    print_stations_rec(T, curr->left);

  printf("%d:\n  ", curr->key);
  int i;
  for (i = 0; i < curr->vehicles->size; i++)
    printf("%d ", curr->vehicles->heap[i]);
  printf("\n");

  if (curr->right != T->nil)
    print_stations_rec(T, curr->right);
}
void print_stations(station_tree_t *T) {
  return print_stations_rec(T, T->root);
}

#pragma endregion

#pragma region commands

void cmd_add_station(station_tree_t *T) {
  int dist, vn;
  scanf("%d %d", &dist, &vn);

  if (find_station(T, dist) != T->nil) {
    printf("non aggiunta\n");
    char *_ = "";
    scanf("%[^\n]", _); // Discards the rest of the line
  } else {
    int i;
    int *vehicles = malloc(sizeof(int) * vn);

    if (vehicles) {
      for (i = 0; i < vn; i++)
        scanf("%d", &vehicles[i]);

      station_t *s = create_station(T, dist, vehicles, vn);

      insert_station(T, s);

      printf("aggiunta\n");
    } else
      printf("[cmd_add_station] Allocation error.\n");
  }
}

#pragma endregion
