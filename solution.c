#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct vehicle_ {
  int range;

  struct vehicle_ *left, *right;
} vehicle_t;

typedef struct station_ {
  int distance;
  vehicle_t *vehicles;

  struct station_ *left, *right, *parent;
  int red;
} station_t;

station_t *cmd_add_station(station_t *stations_root, int distance,
                           int vehicle_n, int vehicles[]);

int parent(int i);
int left(int i);
int right(int i);
vehicle_t *build_vehicle_heap(int vehicles_n, int *ranges);
vehicle_t *build_vehicle_heap_node(int vehicles_n, int *ranges, int curr);
station_t *find_station(station_t *stations_root, int distance);
void heapSort(int *, int);
station_t *insert_station(station_t *root, int distance, int vehicle_n,
                          int *vehicles);
void insert_station_fixup(station_t *root, station_t *z);
void max_heapify(int vehicles_n, int *ranges, int i);
void rotate_left(station_t *root, station_t *x);
void rotate_right(station_t *root, station_t *x);

void swap(int *a, int *b) {

  int temp = *a;
  *a = *b;
  *b = temp;
}

int main() {
  station_t *stations_root = NULL;
  int i;

  char command[20];
  // eof_res should be -1 when there is no stdin left to read
  int eof_res = scanf("%s", command);

  while (eof_res >= 0) {
    if (!strcmp(command, "aggiungi-stazione")) {
      int distance, vehicle_n;
      scanf("%d %d", &distance, &vehicle_n);

      int *ranges = malloc(sizeof(int) * vehicle_n);

      if (ranges) {
        for (i = 0; i < vehicle_n; i++)
          scanf("%d", &ranges[i]);

        stations_root =
            cmd_add_station(stations_root, distance, vehicle_n, ranges);
      } else {
        printf("[main] Allocation error.\n");
        return -1;
      }
    } else if (!strcmp(command, "demolisci-stazione")) {
    } else if (!strcmp(command, "aggiungi-auto")) {
    } else if (!strcmp(command, "rottama-auto")) {
    } else if (!strcmp(command, "pianifica-percorso")) {
    } else {
      printf("Unknown command: %s\n", command);
      return -1;
    }

    eof_res = scanf("%s", command);
  }
}

#pragma region utilities
// I'm converting from the 1-based pseudocode to 0-based C arrays
int parent(int i) { return ((i + 1) / 2) - 1; }
int left(int i) { return (2 * (i + 1)) - 1; }
int right(int i) { return ((2 * (i + 1)) + 1) - 1; }

vehicle_t *build_vehicle_heap(int vehicles_n, int *ranges) {
  int i;

  heapSort(ranges, vehicles_n);

  return build_vehicle_heap_node(vehicles_n, ranges, 0);
}

vehicle_t *build_vehicle_heap_node(int vehicles_n, int *ranges, int curr) {
  int l = left(curr), r = right(curr);

  vehicle_t *node = malloc(sizeof(vehicle_t));
  if (node) {
    node->range = ranges[curr];

    if (l < vehicles_n)
      node->left = build_vehicle_heap_node(vehicles_n, ranges, l);
    else
      node->left = NULL;
    if (r < vehicles_n)
      node->left = build_vehicle_heap_node(vehicles_n, ranges, r);
    else
      node->left = NULL;
  } else
    printf("[build_vehicle_heap_node] Allocation error.");

  return node;
}

station_t *find_root(station_t *node) {
  if (node && node->parent != NULL)
    return find_root(node->parent);
  else
    return node;
}

station_t *find_station(station_t *stations_root, int distance) {
  if (stations_root == NULL)
    return NULL;

  int curr = stations_root->distance;
  if (curr == distance)
    return stations_root;

  if (curr > distance)
    return find_station(stations_root->left, distance);
  else
    return find_station(stations_root->right, distance);
}

void heapSort(int arr[], int N) {
  int i;
  // Build max heap
  for (i = N / 2 - 1; i >= 0; i--)

    max_heapify(N, arr, i);

  // Heap sort
  for (i = N - 1; i >= 0; i--) {

    swap(&arr[0], &arr[i]);

    // Heapify root element
    // to get highest element at
    // root again
    max_heapify(i, arr, 0);
  }
}
station_t *insert_station(station_t *root, int distance, int vehicle_n,
                          int *vehicles) {
  station_t *y = NULL;
  station_t *x = root;

  station_t *z = malloc(sizeof(station_t));
  if (z) {
    z->parent = NULL;
    z->distance = distance;
    z->vehicles = build_vehicle_heap(vehicle_n, vehicles);

    if (root == NULL)
      return z;

    while (x != NULL) {
      y = x;
      if (z->distance < x->distance)
        x = x->left;
      else
        x = x->right;
    }

    z->parent = y;
    if (y == NULL)
      root = z;
    else if (z->distance < y->distance)
      y->left = z;
    else
      y->right = z;

    z->left = NULL;
    z->right = NULL;
    z->red = 1;

    insert_station_fixup(root, z);
  } else
    printf("[insert_station] Allocation error.");

  return find_root(z);
}

void insert_station_fixup(station_t *root, station_t *z) {
  station_t *x, *y;

  if (root == z)
    root->red = 0;
  else {
    x = z->parent;

    if (x->red) {
      if (x == x->parent->left) {
        y = x->parent->right;

        if (y->red) {
          x->red = 0;
          y->red = 0;
          x->parent->red = 1;
          insert_station_fixup(root, x->parent);
        } else if (z == x->right) {
          z = x;
          rotate_left(root, z);
          x = z->parent;
        }

        x->red = 0;
        x->parent->red = 1;
        rotate_right(root, x->parent);
      } else {
        y = x->parent->left;

        if (y->red) {
          x->red = 0;
          y->red = 0;
          x->parent->red = 1;
          insert_station_fixup(root, x->parent);
        } else if (z == x->left) {
          z = x;
          rotate_right(root, z);
          x = z->parent;
        }

        x->red = 0;
        x->parent->red = 1;
        rotate_left(root, x->parent);
      }
    }
  }
}

void max_heapify(int vehicles_n, int *ranges, int i) {
  int l = left(i), r = right(i);
  int max;

  if (l < vehicles_n && ranges[l] > ranges[i])
    max = l;
  else
    max = i;
  if (r < vehicles_n && ranges[r] > ranges[max])
    max = r;

  if (max != i) {
    int tmp = ranges[i];
    ranges[i] = ranges[max];
    ranges[max] = tmp;
    max_heapify(vehicles_n, ranges, max);
  }
}

void rotate_left(station_t *root, station_t *x) {
  station_t *y = x->right;
  x->right = y->left;

  if (y->left != NULL)
    y->left->parent = x;
  y->parent = x->parent;

  if (x->parent == NULL)
    root = y;
  else if (x == x->parent->left)
    x->parent->left = y;
  else
    x->parent->right = y;

  y->left = x;
  x->parent = y;
}

void rotate_right(station_t *root, station_t *x) {
  station_t *y = x->left;
  x->left = y->right;

  if (y->right != NULL)
    y->right->parent = x;
  y->parent = x->parent;

  if (x->parent == NULL)
    root = y;
  else if (x == x->parent->right)
    x->parent->right = y;
  else
    x->parent->left = y;

  y->right = x;
  x->parent = y;
}

// vehicle_t *insert_vehicle(vehicle_t *vehicles_root, int range) {

// }
#pragma endregion

#pragma region commands

station_t *cmd_add_station(station_t *stations_root, int distance,
                           int vehicle_n, int vehicles[]) {
  if (find_station(stations_root, distance) != NULL)
    printf("non aggiunta\n");
  else {
    stations_root =
        insert_station(stations_root, distance, vehicle_n, vehicles);
    printf("aggiunta\n");
  }

  return stations_root;
}

#pragma endregion
