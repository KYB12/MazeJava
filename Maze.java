
import java.util.ArrayList;
import java.util.Collections;
import java.util.ArrayDeque;
import java.util.HashSet;
import javalib.impworld.*;
import java.awt.Color;
import javalib.worldimages.*;
import tester.Tester;
import java.util.HashMap;
import java.util.Random;

// Represents a non border vertex. Edges represents the walls
class Vertex {
  // Coordinates of the vertex
  private final int x;
  private final int y;

  // Mini details
  private final int size;
  private final int wall;

  // Only contains edges that have a direct path to another node; no walls
  final ArrayList<Edge> outEdges;

  Vertex(int x, int y, int wall, int size, ArrayList<Edge> outEdges) {
    this.x = x;
    this.y = y;

    this.size = size;
    this.wall = wall;
    this.outEdges = outEdges;
  }

  Vertex(int x, int y, int wall, int size) {
    this(x, y, wall, size, new ArrayList<>());
  }

  // EFFECT: adds an edge to the given pixel
  void addEdge(Edge edge) {
    this.outEdges.add(edge);
  }

  // returns true if there is only one edge
  boolean isDeadEnd() {
    return this.outEdges.size() == 1;
  }

  // EFFECT: Colors the pixels of this vertex in the given image given color
  void render(ComputedPixelImage pixelImage, Color c) {
    ComputedPixelImage image = pixelImage;
    Color toColor = c;

    int startX = 12 * this.x + 1;
    int startY = 12 * this.y + 1;
    // iterate through the rows of the vertex it needs to color including its walls
    for (int yAdj = 0; yAdj < this.size; yAdj += 1) {
      // iterate through columns of the vertex it needs to color including its walls
      for (int xAdj = 0; xAdj < this.size; xAdj += 1) {
        int newX = startX + xAdj;
        int newY = startY + yAdj;
        image.setColorAt(newX, newY, toColor);
      }
    }

    // iterate through the edges(walls) and replace the walls with a specified color
    for (Edge e : this.outEdges) {
      e.knockDown(image, toColor, this);
    }

  }

  // EFFECT: Erases the wall in a given direction in the given image and replaces
  // with given color
  // Throws error if direction is invalid ("Top", "Bottom", "Left", "Right");
  void erase(String direction, ComputedPixelImage image, Color c) {
    // start coordinates of pixel
    int xCoord = this.x * (size + (2 * this.wall));
    int maxX = xCoord + this.wall + this.size;
    int yCoord = this.y * (size + (2 * this.wall));
    int maxY = yCoord + this.wall + this.size;
    if (direction.equals("Top")) {
      this.eraseHelp(xCoord + this.wall, maxX, yCoord, true, image, c);
    }
    else if (direction.equals("Bottom")) {
      this.eraseHelp(xCoord + this.wall, maxX, maxY, true, image, c);
    }
    else if (direction.equals("Left")) {
      this.eraseHelp(yCoord + this.wall, maxY, xCoord, false, image, c);
    }
    else if (direction.equals("Right")) {
      this.eraseHelp(yCoord + this.wall, maxY, maxX, false, image, c);
    }
    else {
      throw new IllegalArgumentException(direction + " is not a valid direction. Try again with"
          + " \"Top\", \"Bottom\", \"Left\", or \"Right\"");
    }
  }

  // Erases a wall on the given image given a starting x and y and a constant to
  // adjust it by T - horizontal wall F - vertical wall
  void eraseHelp(int min, int max, int constant, boolean mode, ComputedPixelImage image, Color c) {
    // Iterate through the size of the wall
    for (int wallSize = 0; wallSize < this.wall; wallSize += 1) {
      // Iterate through the minimum and maximum points
      for (int place = min; place < max; place += 1) {
        if (mode) {
          image.setPixel(place, constant + wallSize, c);
        }
        else {
          image.setPixel(constant + wallSize, place, c);
        }
      }
    }
  }

  // T - head F - tail
  // Adds the children of this node to the worklist and to path first or last
  // based on mode
  void addTo(boolean mode, ArrayDeque<Vertex> worklist, ArrayList<Vertex> path,
      ArrayDeque<ArrayList<Vertex>> minerPath, HashSet<Vertex> visited) {
    // Go in reverse order instead of for each since this works better for our
    // representation
    for (int edge = this.outEdges.size() - 1; edge >= 0; edge -= 1) {
      this.outEdges.get(edge).addTo(mode, worklist, path, minerPath, this, visited);
    }
  }

  // Determins if there is a direct path to the another vertex
  boolean hasPathTo(Vertex other) {
    // Iterate through the edges of the vertex
    for (Edge e : this.outEdges) {
      // We used an if statement instead of just returning since we wanted the loop to
      // keep going even if first edge returned false
      if (e.containsBoth(this, other)) {
        return true;
      }
    }
    return false;
  }
}

//Represents the connections between nodes. 
class Edge implements Comparable<Edge> {
  private final Vertex from;
  private final Vertex to;
  private final int weight;
  // True = vertical, False = horizontal
  private final boolean direction;

  Edge(Vertex from, Vertex to, int weight, boolean direction) {
    this.from = from;
    this.to = to;
    this.weight = weight;
    this.direction = direction;
  }

  public int compareTo(Edge other) {
    return this.weight - other.weight;
  }

  // EFFECT: Removes the wall between the from and to vertices in the given image
  void knockDown(ComputedPixelImage image, Color c, Vertex v) {
    if (this.from == v) {
      if (direction) {
        v.erase("Bottom", image, c);
      }
      else {
        v.erase("Right", image, c);
      }
    }
    else {
      if (direction) {
        v.erase("Top", image, c);
      }
      else {
        v.erase("Left", image, c);
      }
    }
  }

  void knockDownBoth(ComputedPixelImage theImage, Color c) {
    this.knockDown(theImage, c, from);
    this.knockDown(theImage, c, to);
  }

  // EFFECT: If this edge's vertices have the same representative, then does
  // nothing. Otherwise, sets one representative as the representative of the
  // other vertex's representative. Also removes the edge from the two nodes and
  // adds the edge to the hashset given
  // ACCUMULATOR: Adds 1 to the given count if something was reassigned, and adds
  // 0 otherwise
  int reassign(HashMap<Vertex, Vertex> map, HashSet<Edge> edges, int count) {
    // Gets the representatives of each vertex
    Vertex leftVert = this.from;
    Vertex rightVert = this.to;

    // go through the hashmap until we find a Vertex whose representative is itself
    while (map.get(leftVert) != leftVert) {
      leftVert = map.get(leftVert);
    }
    // go through the hashmap until we find a Vertex whose representative is itself
    while (map.get(rightVert) != rightVert) {
      rightVert = map.get(rightVert);
    }

    // They don't have same rep so it is a wall to be removed
    if (leftVert != rightVert) {
      // unionize
      map.put(leftVert, rightVert);
      edges.add(this);
      this.from.addEdge(this);
      this.to.addEdge(this);
      return count + 1;
    }
    else {
      return count;
    }
  }

  // T - head(DFS) F- tail(BFS)
  // Adds the vertices in this edge to the worklist and path except if its the
  // same as the given one
  void addTo(boolean mode, ArrayDeque<Vertex> worklist, ArrayList<Vertex> path,
      ArrayDeque<ArrayList<Vertex>> minerPath, Vertex duplicate, HashSet<Vertex> visited) {
    // Make shallow copies of the path to have 2 seperate paths; one with from
    // vertex and one with to vertex

    ArrayList<Vertex> dummy1 = new ArrayList<>();
    // iterate through all the vertices in the path and add them to the first list
    for (Vertex v : path) {
      dummy1.add(v);
    }

    ArrayList<Vertex> dummy2 = new ArrayList<>();
    // iterate through all the vertices in the path and add them to the second list
    for (Vertex v : path) {
      dummy2.add(v);
    }

    dummy1.add(from);
    dummy2.add(to);

    if (mode) {
      if (this.from != duplicate && !visited.contains(from)) {
        worklist.addFirst(from);
        minerPath.addFirst(dummy1);
      }
      if (this.to != duplicate && !visited.contains(to)) {
        worklist.addFirst(to);
        minerPath.addFirst(dummy2);
      }
    }
    else {
      if (this.from != duplicate && !visited.contains(from)) {
        worklist.addLast(from);
        minerPath.addLast(dummy1);
      }
      if (this.to != duplicate && !visited.contains(to)) {
        worklist.addLast(to);
        minerPath.addLast(dummy2);
      }
    }
  }

  // Determines if both vertices are in the edge regardless of order
  boolean containsBoth(Vertex first, Vertex second) {
    return (this.from == first && this.to == second) || (this.from == second && this.to == first);
  }

}

class Maze {
  private final int width;
  private final int height;

  // Not final since the start and end can be changed based on preference for our
  // implementation of maze
  private Vertex start;
  private Vertex end;

  // List of all the vertices
  private final ArrayList<Vertex> vertices;
  // List of all the visited vertices
  private final HashSet<Vertex> visited;

  // The list of edges making up the MST
  private final HashSet<Edge> tree;
  // Represents each vertex and its representative
  private final HashMap<Vertex, Vertex> unionFind;

  // Vertices to analyze
  private final ArrayDeque<Vertex> worklist;
  // Paths to analyze
  private final ArrayDeque<ArrayList<Vertex>> minerPath;
  // Holds the vertices that make up the solution
  private final ArrayList<Vertex> solutionPath;

  // Heap of all the edges
  private final ArrayList<Edge> edgesArray;

  // Not final since the distance from the origin changes
  private int distance;

  // Not final since the miner moves around. Vertex and Maze can edit it since
  // they know whether or not a miner can move or not
  private int minerLocation;

  // Image representation of the Maze
  private final ComputedPixelImage theImage;
  // Heatmap of distance from start
  private final ComputedPixelImage redBlue;
  // Heatmap of distance from end
  private final ComputedPixelImage blueRed;

  // Different colors for different states
  private final Color frontier;
  private final Color vis;
  private final Color current;
  private final Color solution;

  private int largest;

  // Modes:
  // 0 - Normal
  // 1 - horizontal corridor bias
  // 2 - vertical corridor bias
  Maze(int width, int height, int wall, int size, ArrayList<Edge> edges, ArrayList<Vertex> vertices,
      ArrayList<Vertex> solution, int mode) {
    if (width < 0 || height < 0) {
      throw new IllegalArgumentException("Size must be greater than 0");
    }
    this.distance = 0;
    this.largest = 0;
    ArrayList<Vertex> dummyList = new ArrayList<>();
    ArrayList<Edge> dummyHeap = new ArrayList<>();
    this.tree = new HashSet<>();
    this.unionFind = new HashMap<>();
    this.visited = new HashSet<>();
    this.minerPath = new ArrayDeque<>();
    this.solutionPath = solution;
    this.width = width;
    this.height = height;
    this.minerLocation = 0;
    int computedSize = ((2 * wall) + size);
    this.theImage = new ComputedPixelImage(width * computedSize, height * computedSize);
    this.redBlue = new ComputedPixelImage(width * computedSize, height * computedSize);
    this.blueRed = new ComputedPixelImage(width * computedSize, height * computedSize);
    this.frontier = Color.lightGray;
    this.vis = Color.CYAN;
    this.current = Color.BLUE;
    this.solution = Color.GREEN;

    // iterate through every vertex that needs to be created
    for (int vertex = 0; vertex < width * height; vertex += 1) {
      // random number between 0 - 10
      int randEdge = (int) (new Random().nextDouble() * 10);
      int x = vertex % width;
      int y = vertex / width;
      Vertex current = new Vertex(x, y, wall, size);

      // not on topRow:
      if (vertex >= width) {
        if (mode == 1) {
          // give vertical walls lower weights so they get removed first
          randEdge = randEdge - 100;
        }
        Edge toAdd = new Edge(dummyList.get(vertex - width), current, randEdge, true);
        dummyHeap.add(toAdd);
      }

      // not on left edge:
      randEdge = (int) (new Random().nextDouble() * 10);
      if (vertex % width > 0) {
        if (mode == 2) {
          // give vertical walls lower weights so they get removed first
          randEdge = randEdge - 100;
        }
        Edge toAdd = new Edge(dummyList.get(vertex - 1), current, randEdge, false);
        dummyHeap.add(toAdd);
      }

      int pixelX = computedSize * x;
      int pixelY = computedSize * y;

      current.render(this.theImage, this.frontier);

      // Set up the walls of the pixel
      // Iterate through the size of the pixel horizontally
      for (int colAdjust = 0; colAdjust < computedSize; colAdjust += 1) {
        this.theImage.setColorAt(pixelX + colAdjust, pixelY, Color.darkGray);
        this.theImage.setColorAt(pixelX + colAdjust, pixelY + computedSize - 1, Color.darkGray);

        // For redBlue
        this.redBlue.setColorAt(pixelX + colAdjust, pixelY, Color.darkGray);
        this.redBlue.setColorAt(pixelX + colAdjust, pixelY + computedSize - 1, Color.darkGray);

        // For blueRed
        this.blueRed.setColorAt(pixelX + colAdjust, pixelY, Color.darkGray);
        this.blueRed.setColorAt(pixelX + colAdjust, pixelY + computedSize - 1, Color.darkGray);
      }

      // Iterate through the size of the pixel vertically
      for (int rowAdjust = 0; rowAdjust < computedSize; rowAdjust += 1) {
        this.theImage.setColorAt(pixelX, pixelY + rowAdjust, Color.darkGray);
        this.theImage.setColorAt(pixelX + computedSize - 1, pixelY + rowAdjust, Color.darkGray);

        this.redBlue.setColorAt(pixelX, pixelY + rowAdjust, Color.darkGray);
        this.redBlue.setColorAt(pixelX + computedSize - 1, pixelY + rowAdjust, Color.darkGray);

        this.blueRed.setColorAt(pixelX, pixelY + rowAdjust, Color.darkGray);
        this.blueRed.setColorAt(pixelX + computedSize - 1, pixelY + rowAdjust, Color.darkGray);
      }

      dummyList.add(current);

      // Testing purposes
      if (vertices.size() > 0) {
        this.unionFind.put(vertices.get(vertex), vertices.get(vertex));
      }
      else {
        this.unionFind.put(current, current);
      }

    }

    // Testing purposes
    if (vertices.size() > 0) {
      this.vertices = vertices;
    }
    else {
      this.vertices = dummyList;
    }
    this.start = this.vertices.get(0);
    this.end = this.vertices.get(height * width - 1);
    this.worklist = new ArrayDeque<Vertex>();
    this.worklist.add(this.start);

    // testing purposes
    if (edges.size() > 0) {
      this.edgesArray = edges;
    }
    else {
      this.edgesArray = dummyHeap;
    }
    ArrayList<Vertex> path = new ArrayList<>();
    path.add(start);
    this.minerPath.add(path);
  }

  Maze(int width, int height, int wall, int size, int mode) {
    this(width, height, wall, size, new ArrayList<Edge>(), new ArrayList<Vertex>(),
        new ArrayList<Vertex>(), mode);
  }

  // returns the image of the current maze
  // 0 - BFS
  // 1 - DFS
  // 2- manual
  // 3 - solution
  // 4 - heatmap from start
  // 5 - heatmap from end
  ComputedPixelImage render(int mode) {
    // Iterate through the visited vertecies
    for (Vertex v : this.visited) {
      v.render(theImage, this.vis);
    }

    if ((this.vertices.get(minerLocation) == this.end) || (mode == 3)) {
      // if its BFS, the worklist is the heads and thus are visited
      // if its DFS, only the head is visited and at this point, the maze is solved so
      // there is nothing to color
      if (mode == 0) {
        for (Vertex v : this.worklist) {
          v.render(theImage, this.vis);
        }
      }
      // Itetate through the solution vertices
      for (Vertex v : this.solutionPath) {
        v.render(theImage, this.solution);
      }
    }
    else {
      if (mode == 0) {
        // Go through the vorklist since breadth first does one layer at a time
        for (Vertex v : this.worklist) {
          v.render(theImage, this.current);
        }
      }
      else if (mode == 1) {
        // Only get the head of the worklist since dfs does one vertex at a time
        if (!this.worklist.isEmpty()) {
          this.worklist.getFirst().render(theImage, this.current);
        }
      }
      else if (mode == 2) {
        this.vertices.get(this.minerLocation).render(theImage, this.current);
      }
      else if (mode == 4) {
        return this.redBlue;
      }
      else if (mode == 5) {
        return this.blueRed;
      }
    }

    return this.theImage;
  }

  // Constructs the maze
  void construct() {
    int edgeCount = 0;
    int numVert = this.vertices.size();
    // iterates through heap until we have the amount of edges needed
    Collections.sort(edgesArray);
    while (edgeCount < numVert - 2) {
      Edge toRemove = this.edgesArray.remove(0);
      edgeCount = toRemove.reassign(this.unionFind, this.tree, edgeCount);
    }

    // For every path, knock down the wall
    for (Edge e : this.tree) {
      e.knockDownBoth(this.theImage, Color.lightGray);
      // same for redblue
      e.knockDownBoth(this.redBlue, Color.lightGray);
      // Same for bluered
      e.knockDownBoth(this.blueRed, Color.lightGray);
    }

  }

  // T - DFS F - BFS
  // Traverse the maze DFS or BFS and determines whether or not the end of the
  // maze has been reached and finds that path
  // EFFECT: Alters the worklist and minerpath to remove and add to them
  // accordingly to BFS or DFS order and adds the solution to solutionpath
  boolean traverse(boolean mode) {
    if (!mode) {
      int max = this.worklist.size();
      // the initial worklist is all the vertices that need to be processed. So
      // iterate that many times to go one layer at a time
      for (int i = 0; i < max; i += 1) {
        if (!this.worklist.isEmpty() && this.traverseHelp(mode)) {
          return true;
        }
      }
      return false;
    }
    else {
      // Iterates the worklist while the head is not a dead end or start or end.
      // This way, we keep going down one "depth" at a time rather than one node at a
      // time
      while (!this.worklist.getFirst().isDeadEnd() || this.worklist.getFirst() == this.start
          || this.worklist.getFirst() == this.end) {
        if (this.traverseHelp(mode)) {
          return true;
        }
      }
      // If we reach a dead end, remove it from the path and worklist
      this.visited.add(this.worklist.remove());
      this.minerPath.remove();
      return false;
    }
  }

  // Returns whether or not the bottom right vertex has been reached
  // EFFECT: Removes the head of the path and worklsit and adds its neighbors to
  // either the head or tail of the minerPath and Worklist based on mode
  // T - DFS F - BFS
  boolean traverseHelp(boolean mode) {
    Vertex current = worklist.remove();
    ArrayList<Vertex> path = this.minerPath.remove();
    if (current == this.end) {
      this.minerLocation = this.vertices.size() - 1;
      current.addTo(mode, this.worklist, path, this.minerPath, this.visited);
      // Add everything from the path to the solution path since it is the solution
      for (Vertex v : path) {
        this.solutionPath.add(v);
      }
      return true;
    }
    if (!this.visited.contains(current)) {
      this.visited.add(current);
      current.addTo(mode, this.worklist, path, this.minerPath, this.visited);
    }
    return false;
  }

  // EFFECT: updates the largest distance to the path with the largest distance
  // Also colors the vertex based on how far it is from the starting point if the
  // boolean is true
  void heatMap(ComputedPixelImage img, boolean color) {
    int totalNum = this.vertices.size();
    // While every node has hasn't been visited
    while (this.visited.size() < totalNum - 1) {
      // Estimate what the maximum path of the
      if (color) {
        float blue = (float) (this.distance) / (float) (this.largest);
        // just in case it goes over
        if (blue > 1) {
          blue = 1;
        }
        Color toColor = new Color(1 - blue, 0, blue);
        // Paint the vertex
        for (Vertex v : this.worklist) {
          v.render(img, toColor);
        }
      }
      this.traverse(false);
      distance += 1;
    }

    // Reset mutated lists
    this.worklist.clear();
    this.worklist.add(start);
    this.minerPath.clear();
    ArrayList<Vertex> path = new ArrayList<>();
    path.add(start);
    minerPath.add(path);
    this.visited.clear();
    this.minerLocation = 0;
    this.largest = this.distance;
    this.distance = 0;
  }

  // Effect: Colors the redBlue image based on how far it from the start. Red
  // means close and blue means far away
  void redBlue() {
    this.colorHeat(this.start, this.redBlue);
  }

  // Effect: Colors the blueRed image based on how far it is from the end. Red
  // means close and blue means far away
  void blueRed() {
    this.colorHeat(this.end, this.blueRed);
  }

  // Effect: Given a starting node, colors every vertex in the image based on how
  // far away it is from that node. Red means close and blue means far away
  void colorHeat(Vertex v, ComputedPixelImage i) {
    this.worklist.clear();
    this.worklist.add(v);
    this.distance = 0;
    this.largest = 0;
    this.heatMap(i, false);
    this.worklist.clear();
    this.worklist.add(v);
    this.heatMap(i, true);
  }

  // finds the solution of this maze
  void findSolution() {
    // traverse maze DFS until a solution is found
    while (!this.traverse(false)) {
      // Do nothing since traverse already finds the path
    }
    // reset the deques and arraylists to original state
    this.minerPath.clear();
    this.visited.clear();
    this.worklist.clear();
    ArrayList<Vertex> path = new ArrayList<>();
    path.add(this.start);
    this.minerPath.add(path);
    this.worklist.add(start);
    this.minerLocation = 0;
  }

  // Move in a specified direction if possible
  // EFFECT: Clears the worklsit as manual doesn't allow automated traversal
  // anymore also changes the minerlocation to the new location.
  public void move(String direction) {
    if (direction.equals("right") && this.minerLocation % this.width < this.width - 1) {
      int newPlace = this.minerLocation + 1;
      this.moveIfCan(newPlace);
    }

    else if (direction.equals("left") && this.minerLocation % this.width > 0) {
      int newPlace = this.minerLocation - 1;
      this.moveIfCan(newPlace);
    }

    else if (direction.equals("up") && this.minerLocation / this.height > 0) {
      int newPlace = this.minerLocation - this.width;
      this.moveIfCan(newPlace);

    }
    else if (direction.equals("down") && this.minerLocation / this.height < this.height - 1) {
      int newPlace = this.minerLocation + this.width;
      this.moveIfCan(newPlace);
    }
    else {
      // do nothing
      return;
    }
  }

  // Effect: Moves the miner to the given location if it can move that way. i.e.
  // no walls blocking it
  void moveIfCan(int location) {
    this.worklist.clear();
    Vertex current = this.vertices.get(this.minerLocation);
    Vertex next = this.vertices.get(location);
    if (current.hasPathTo(next)) {
      this.visited.add(current);
      this.minerLocation = location;
    }
    else {
      // do nothing
      return;
    }
  }
}

//Reprsents the interface that users interact with the maze
class MyWorld extends World {
  // can't make final because the maze changes based on input. The heap gets used
  // up so we have to reinitialize a new maze
  Maze maze;

  // Not final since the user can choose between breadth first or depth first
  // 0 = breadth, 1 = depth, 2 = manual
  private int method;

  // Whether or not to traverse
  private boolean traverse;

  private final int width;
  private final int height;

  // determines whether or not manual input is allowed
  private boolean prohibit;
  // same thing for automated search
  private boolean prohibitDig;

  MyWorld(Maze maze, int width, int height) {
    this.maze = maze;
    this.width = width;
    this.height = height;
    this.traverse = false;
    this.prohibit = false;
    this.prohibitDig = false;
    this.method = 2;
    this.maze.construct();
//    this.maze.findSolution();
//    this.maze.redBlue();
//    this.maze.blueRed();
  }

  // draw out the image scene
  public WorldScene makeScene() {
    WorldScene canvas = new WorldScene(1200, 720);
    canvas.placeImageXY(maze.render(method), 600, 360);
    return canvas;
  }

  // Make a new maze
  public void onKeyEvent(String key) {
    if (key.equals("x")) {
      this.method = 2;
      this.maze = new Maze(width, height, 1, 10, 0);
      this.maze.construct();
      this.maze.findSolution();
      this.maze.redBlue();
      this.maze.blueRed();
      this.traverse = false;
      this.prohibit = false;
      this.prohibitDig = false;
    }
    else if (key.equals("up") && !prohibit) {
      this.method = 2;
      this.prohibitDig = true;
      this.maze.move(key);
    }
    else if (key.equals("down") && !prohibit) {
      this.method = 2;
      this.prohibitDig = true;
      this.maze.move(key);
    }
    else if (key.equals("left") && !prohibit) {
      this.method = 2;
      this.prohibitDig = true;
      this.maze.move(key);
    }
    else if (key.equals("right") && !prohibit) {
      this.method = 2;
      this.prohibitDig = true;
      this.maze.move(key);
    }
    else if (key.equals("b") && !prohibitDig) {
      this.method = 0;
      this.traverse = true;
      this.prohibit = true;
    }
    else if (key.equals("d") && !prohibitDig) {
      this.method = 1;
      this.traverse = true;
      this.prohibit = true;
    }
//    else if (key.equals("s")) {
//      traverse = false;
//      this.method = 4;
//    }
//    else if (key.equals("e")) {
//      traverse = false;
//      this.method = 5;
//    }
    // Horizontal corridor bias
    else if (key.equals("h")) {
      this.method = 2;
      this.maze = new Maze(width, height, 1, 10, 1);
      this.maze.construct();
      this.maze.findSolution();
      this.maze.redBlue();
      this.maze.blueRed();
      this.traverse = false;
      this.prohibit = false;
      this.prohibitDig = false;
    }
    // Vertical corridor bias
    else if (key.equals("v")) {
      this.method = 2;
      this.maze = new Maze(width, height, 1, 10, 2);
      this.maze.construct();
      this.maze.findSolution();
      this.maze.redBlue();
      this.maze.blueRed();
      this.traverse = false;
      this.prohibit = false;
      this.prohibitDig = false;
    }
  }

  // Updates the screen every tick based on conditions
  public void onTick() {
    if (traverse) {
      boolean end = false;
      if (this.method == 0) {
        end = this.maze.traverse(false);
      }
      else if (this.method == 1) {
        end = this.maze.traverse(true);

      }
      if (end) {
        this.traverse = false;
      }
    }
  }

}

class ExamplesMaze {

  Color frontier;
  Color visited;
  Color current;
  Color solution;

  Maze maze1;
  Maze examplesMaze1;

  Vertex A1;
  Vertex A2;
  Vertex A3;
  Vertex A4;
  Vertex B1;
  Vertex B2;
  Vertex B3;
  Vertex B4;
  Vertex C1;
  Vertex C2;
  Vertex C3;
  Vertex C4;
  Vertex D1;
  Vertex D2;
  Vertex D3;
  Vertex D4;

  Edge A1A2;
  Edge A2A3;
  Edge A3A4;
  Edge B1B2;
  Edge B2B3;
  Edge B3B4;
  Edge C1C2;
  Edge C2C3;
  Edge C3C4;
  Edge D1D2;
  Edge D2D3;
  Edge D3D4;

  Edge A1B1;
  Edge B1C1;
  Edge C1D1;

  Edge A2B2;
  Edge B2C2;
  Edge C2D2;

  Edge A3B3;
  Edge B3C3;
  Edge C3D3;

  Edge A4B4;
  Edge B4C4;
  Edge C4D4;

  ArrayList<Edge> exampleHeap;

  ComputedPixelImage image1;
  ComputedPixelImage canvas;
  ComputedPixelImage emptyGrid;
  ComputedPixelImage knockedDownMaze;

  // List of vertices A1 - D4
  ArrayList<Vertex> dummyList;
  // Heap of edges in solution
  ArrayList<Edge> dummyHeap;

  HashSet<Edge> edges;

  ArrayDeque<Vertex> actualWorklist;
  ArrayList<Vertex> actualPath;
  ArrayList<Vertex> solutionPath;
  ArrayDeque<ArrayList<Vertex>> actualMinerPath;

  ArrayDeque<Vertex> expectedWorklist;
  ArrayList<Vertex> expectedPath;
  ArrayDeque<ArrayList<Vertex>> expectedMinerPath;

  HashMap<Vertex, Vertex> actualMap;
  HashMap<Vertex, Vertex> expectedMap;

  IllegalArgumentException direction;

  void init() {
    this.maze1 = new Maze(100, 60, 1, 10, 0);
    this.image1 = new ComputedPixelImage(48, 48);
    this.canvas = new ComputedPixelImage(48, 48);

    A1 = new Vertex(0, 0, 1, 10);
    A2 = new Vertex(1, 0, 1, 10);
    A3 = new Vertex(2, 0, 1, 10);
    A4 = new Vertex(3, 0, 1, 10);
    B1 = new Vertex(0, 1, 1, 10);
    B2 = new Vertex(1, 1, 1, 10);
    B3 = new Vertex(2, 1, 1, 10);
    B4 = new Vertex(3, 1, 1, 10);
    C1 = new Vertex(0, 2, 1, 10);
    C2 = new Vertex(1, 2, 1, 10);
    C3 = new Vertex(2, 2, 1, 10);
    C4 = new Vertex(3, 2, 1, 10);
    D1 = new Vertex(0, 3, 1, 10);
    D2 = new Vertex(1, 3, 1, 10);
    D3 = new Vertex(2, 3, 1, 10);
    D4 = new Vertex(3, 3, 1, 10);

    dummyList = new ArrayList<>();
    dummyList.add(A1);
    dummyList.add(A2);
    dummyList.add(A3);
    dummyList.add(A4);
    dummyList.add(B1);
    dummyList.add(B2);
    dummyList.add(B3);
    dummyList.add(B4);
    dummyList.add(C1);
    dummyList.add(C2);
    dummyList.add(C3);
    dummyList.add(C4);
    dummyList.add(D1);
    dummyList.add(D2);
    dummyList.add(D3);
    dummyList.add(D4);

    A1A2 = new Edge(A1, A2, 1, false);
    A2A3 = new Edge(A2, A3, 2, false);
    A3A4 = new Edge(A3, A4, 1, false);
    B1B2 = new Edge(B1, B2, 1, false);
    B2B3 = new Edge(B2, B3, 1, false);
    B3B4 = new Edge(B3, B4, 1, false);
    C1C2 = new Edge(C1, C2, 1, false);
    C2C3 = new Edge(C2, C3, 1, false);
    C3C4 = new Edge(C3, C4, 1, false);
    D1D2 = new Edge(D1, D2, 1, false);
    D2D3 = new Edge(D2, D3, 1, false);
    D3D4 = new Edge(D3, D4, 1, false);

    A1B1 = new Edge(A1, B1, 1, true);
    B1C1 = new Edge(B1, C1, 1, true);
    C1D1 = new Edge(C1, D1, 1, true);

    A2B2 = new Edge(A2, B2, 1, true);
    B2C2 = new Edge(B2, C2, 1, true);
    C2D2 = new Edge(C2, D2, 1, true);

    A3B3 = new Edge(A3, B3, 1, true);
    B3C3 = new Edge(B3, C3, 1, true);
    C3D3 = new Edge(C3, D3, 1, true);

    A4B4 = new Edge(A4, B4, 1, true);
    B4C4 = new Edge(B4, C4, 1, true);
    C4D4 = new Edge(C4, D4, 1, true);

    actualMap = new HashMap<>();
    actualMap.put(A1, A1);
    actualMap.put(A2, A2);
    actualMap.put(A3, A3);
    actualMap.put(A4, A4);
    actualMap.put(B1, B1);
    actualMap.put(B2, B2);
    actualMap.put(B3, B3);
    actualMap.put(B4, B4);
    actualMap.put(C1, C1);
    actualMap.put(C2, C2);
    actualMap.put(C3, C3);
    actualMap.put(C4, C4);
    actualMap.put(D1, D1);
    actualMap.put(D2, D2);
    actualMap.put(D3, D3);
    actualMap.put(D4, D4);

    expectedMap = new HashMap<>();
    expectedMap.put(A1, A1);
    expectedMap.put(A2, A2);
    expectedMap.put(A3, A3);
    expectedMap.put(A4, A4);
    expectedMap.put(B1, B1);
    expectedMap.put(B2, B2);
    expectedMap.put(B3, B3);
    expectedMap.put(B4, B4);
    expectedMap.put(C1, C1);
    expectedMap.put(C2, C2);
    expectedMap.put(C3, C3);
    expectedMap.put(C4, C4);
    expectedMap.put(D1, D1);
    expectedMap.put(D2, D2);
    expectedMap.put(D3, D3);
    expectedMap.put(D4, D4);

    edges = new HashSet<>();
    edges.add(A1A2);
    edges.add(A3A4);
    edges.add(A2B2);
    edges.add(A4B4);
    edges.add(B1B2);
    edges.add(B2B3);
    edges.add(B3B4);
    edges.add(B1C1);
    edges.add(B3C3);
    edges.add(B4C4);
    edges.add(C1D1);
    edges.add(C2D2);
    edges.add(C3D3);
    edges.add(D1D2);
    edges.add(D3D4);

    dummyHeap = new ArrayList<>();
    dummyHeap.add(A1A2);
    dummyHeap.add(A3A4);
    dummyHeap.add(A2B2);
    dummyHeap.add(A4B4);
    dummyHeap.add(B1B2);
    dummyHeap.add(B2B3);
    dummyHeap.add(B3B4);
    dummyHeap.add(B1C1);
    dummyHeap.add(B3C3);
    dummyHeap.add(B4C4);
    dummyHeap.add(C1D1);
    dummyHeap.add(C2D2);
    dummyHeap.add(C3D3);
    dummyHeap.add(D1D2);
    dummyHeap.add(D3D4);
    this.examplesMaze1 = new Maze(4, 4, 1, 10, dummyHeap, dummyList, new ArrayList<>(), 0);

    exampleHeap = new ArrayList<>();
    exampleHeap.add(A1A2);
    exampleHeap.add(A3A4);
    exampleHeap.add(A2B2);
    exampleHeap.add(A4B4);
    exampleHeap.add(B1B2);
    exampleHeap.add(B2B3);
    exampleHeap.add(B3B4);
    exampleHeap.add(B1C1);
    exampleHeap.add(B3C3);
    exampleHeap.add(B4C4);
    exampleHeap.add(C1D1);
    exampleHeap.add(C2D2);
    exampleHeap.add(C3D3);
    exampleHeap.add(D1D2);
    exampleHeap.add(D3D4);

    knockedDownMaze = new Maze(4, 4, 1, 10, dummyHeap, dummyList, new ArrayList<>(), 0).render(2);
    // take down the initial walls while there are edges left
    while (this.exampleHeap.size() > 0) {
      exampleHeap.remove(0).knockDownBoth(knockedDownMaze, Color.lightGray);
    }
    // make the beggining blue
    for (int row = 1; row < 11; row += 1) {
      knockedDownMaze.setColorAt(11, row, Color.blue);
    }

    // make an empty grid
    emptyGrid = new ComputedPixelImage(48, 48);
    emptyGrid = new Maze(4, 4, 1, 10, dummyHeap, dummyList, new ArrayList<>(), 0).render(2);

    this.actualWorklist = new ArrayDeque<>();
    this.actualPath = new ArrayList<>();
    this.actualMinerPath = new ArrayDeque<>();

    this.expectedWorklist = new ArrayDeque<>();
    this.expectedPath = new ArrayList<>();
    this.expectedMinerPath = new ArrayDeque<>();

    this.solutionPath = new ArrayList<>();
    solutionPath.add(A1);
    solutionPath.add(A2);
    solutionPath.add(B2);
    solutionPath.add(B3);
    solutionPath.add(C3);
    solutionPath.add(D3);
    solutionPath.add(D4);

    this.frontier = Color.lightGray;
    this.visited = Color.CYAN;
    this.current = Color.BLUE;
    this.solution = Color.GREEN;

    this.direction = new IllegalArgumentException("Misinput is not a valid direction."
        + " Try again with \"Top\", \"Bottom\", \"Left\", or \"Right\"");
  }

  void testBigBang(Tester t) {
    this.init();
    maze1.construct();
    World myWorld = new MyWorld(examplesMaze1, 100, 60);
    myWorld.bigBang(1200, 720, 0.001);
  }

}
