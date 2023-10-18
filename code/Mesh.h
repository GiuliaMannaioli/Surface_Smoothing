#ifndef MESH_H
#define MESH_H

#include <glad/glad.h>
#include <vector>
#include <memory>
#include <glm/glm.hpp>
#include <glm/ext.hpp>
#include <glm/gtx/vector_angle.hpp>

#include <map>
#include <set>
#include <string>
#include <iostream>


class Mesh {
public:
    virtual ~Mesh();

    const std::vector<glm::vec3>& vertexPositions() const { return _vertexPositions; }
    std::vector<glm::vec3>& vertexPositions() { return _vertexPositions; }

    const std::vector<glm::vec3>& vertexNormals() const { return _vertexNormals; }
    std::vector<glm::vec3>& vertexNormals() { return _vertexNormals; }

    const std::vector<glm::vec2>& vertexTexCoords() const { return _vertexTexCoords; }
    std::vector<glm::vec2>& vertexTexCoords() { return _vertexTexCoords; }

    const std::vector<glm::uvec3>& triangleIndices() const { return _triangleIndices; }
    std::vector<glm::uvec3>& triangleIndices() { return _triangleIndices; }

    /// Compute the parameters of a sphere which bounds the mesh
    void computeBoundingSphere(glm::vec3& center, float& radius) const;

    void recomputePerVertexNormals(bool angleBased = false);
    void recomputePerVertexTextureCoordinates();

    void init();
    void initOldGL();
    void render();
    void clear();

    void addPlan(float square_half_side = 1.0f);

    void subdivideLinear() {
        std::vector<glm::vec3> newVertices = _vertexPositions;
        std::vector<glm::uvec3> newTriangles;

        struct Edge {
            unsigned int a, b;
            Edge(unsigned int c, unsigned int d) : a(std::min<unsigned int>(c, d)), b(std::max<unsigned int>(c, d)) {}
            bool operator < (Edge const& o) const { return a < o.a || (a == o.a && b < o.b); }
            bool operator == (Edge const& o) const { return a == o.a && b == o.b; }
        };
        std::map< Edge, unsigned int > newVertexOnEdge;
        for (unsigned int tIt = 0; tIt < _triangleIndices.size(); ++tIt) {
            unsigned int a = _triangleIndices[tIt][0];
            unsigned int b = _triangleIndices[tIt][1];
            unsigned int c = _triangleIndices[tIt][2];


            Edge Eab(a, b);
            unsigned int oddVertexOnEdgeEab = 0;
            if (newVertexOnEdge.find(Eab) == newVertexOnEdge.end()) {
                newVertices.push_back((_vertexPositions[a] + _vertexPositions[b]) / 2.f);
                //cerca il nuovo vertice al centro del nuovo edge, da usare nella prossima iterazione
                oddVertexOnEdgeEab = newVertices.size() - 1;
                newVertexOnEdge[Eab] = oddVertexOnEdgeEab;
            }
            else { oddVertexOnEdgeEab = newVertexOnEdge[Eab]; }
            //quindi se già lavevo fatto xk è un edge in comune allora tengo quello


            Edge Ebc(b, c);
            unsigned int oddVertexOnEdgeEbc = 0;
            if (newVertexOnEdge.find(Ebc) == newVertexOnEdge.end()) {
                newVertices.push_back((_vertexPositions[b] + _vertexPositions[c]) / 2.f);
                oddVertexOnEdgeEbc = newVertices.size() - 1;
                newVertexOnEdge[Ebc] = oddVertexOnEdgeEbc;
            }
            else { oddVertexOnEdgeEbc = newVertexOnEdge[Ebc]; }

            Edge Eca(c, a);
            unsigned int oddVertexOnEdgeEca = 0;
            if (newVertexOnEdge.find(Eca) == newVertexOnEdge.end()) {
                newVertices.push_back((_vertexPositions[c] + _vertexPositions[a]) / 2.f);
                oddVertexOnEdgeEca = newVertices.size() - 1;
                newVertexOnEdge[Eca] = oddVertexOnEdgeEca;
            }
            else { oddVertexOnEdgeEca = newVertexOnEdge[Eca]; }

            // set new triangles :
            newTriangles.push_back(glm::uvec3(a, oddVertexOnEdgeEab, oddVertexOnEdgeEca));
            newTriangles.push_back(glm::uvec3(oddVertexOnEdgeEab, b, oddVertexOnEdgeEbc));
            newTriangles.push_back(glm::uvec3(oddVertexOnEdgeEca, oddVertexOnEdgeEbc, c));
            newTriangles.push_back(glm::uvec3(oddVertexOnEdgeEab, oddVertexOnEdgeEbc, oddVertexOnEdgeEca));
        }

        // after that:
        _triangleIndices = newTriangles;
        _vertexPositions = newVertices;
        recomputePerVertexNormals();
        recomputePerVertexTextureCoordinates();
    }

    void subdivideLoop() {
        std::vector<glm::vec3> newVertices = _vertexPositions;
        std::vector<glm::uvec3> newTriangles;

        // Here we create the struct that defines an edge
        struct Edge {
            unsigned int a, b;
            Edge(unsigned int c, unsigned int d) : a(std::min<unsigned int>(c, d)), b(std::max<unsigned int>(c, d)) {}
            bool operator < (Edge const& o) const { return a < o.a || (a == o.a && b < o.b); }
            bool operator == (Edge const& o) const { return a == o.a && b == o.b; }
        };

        std::map< Edge, unsigned int > newVertexOnEdge; // this will be useful to find out whether we already inserted an odd vertex or not
        std::vector< std::map<Edge, unsigned int>> trianglesOnEdge(_vertexPositions.size()); ; // this will be useful to find out if an edge is boundary or not
        std::vector< std::set< unsigned int > > neighboringVertices(_vertexPositions.size()); // this will be used to store the adjacent vertices, i.e., neighboringVertices[i]
        // will be the list of vertices that are adjacent to vertex i.

        // I) First, compute the valences of the even vertices, the neighboring vertices 
                        //required to update the position of the even vertices, and the boundaries:

        for (unsigned int tIt = 0; tIt < _triangleIndices.size(); ++tIt) {
            unsigned int a = _triangleIndices[tIt][0];
            unsigned int b = _triangleIndices[tIt][1];
            unsigned int c = _triangleIndices[tIt][2];

            //TODO: Remember the faces shared by the edge
            Edge Eab(a, b);
            Edge Ebc(b, c);
            Edge Eac(a, c);

            //add triangles to the map

              // a
            trianglesOnEdge[a][Eab] = (trianglesOnEdge[a].find(Eab) == trianglesOnEdge[a].end()) ? 1 : 2;
            trianglesOnEdge[a][Eac] = (trianglesOnEdge[a].find(Eac) == trianglesOnEdge[a].end()) ? 1 : 2;
            // b
            trianglesOnEdge[b][Eab] = (trianglesOnEdge[b].find(Eab) == trianglesOnEdge[b].end()) ? 1 : 2;
            trianglesOnEdge[b][Ebc] = (trianglesOnEdge[b].find(Ebc) == trianglesOnEdge[b].end()) ? 1 : 2;
            // c
            trianglesOnEdge[c][Eac] = (trianglesOnEdge[c].find(Eac) == trianglesOnEdge[c].end()) ? 1 : 2;
            trianglesOnEdge[c][Ebc] = (trianglesOnEdge[c].find(Ebc) == trianglesOnEdge[c].end()) ? 1 : 2;



            //to store the neighbors

            neighboringVertices[a].insert(b); //quindi sul vicino di a sono aggiunti b e c
            neighboringVertices[a].insert(c);
            neighboringVertices[b].insert(a);
            neighboringVertices[b].insert(c);
            neighboringVertices[c].insert(a);
            neighboringVertices[c].insert(b);


        }

        // The valence of a vertex is the number of adjacent vertices:
        std::vector< unsigned int > evenVertexValence(_vertexPositions.size(), 0);
        for (unsigned int v = 0; v < _vertexPositions.size(); ++v) {
            evenVertexValence[v] = neighboringVertices[v].size();

        }

        //TODO: Identify even vertices (clue: check the number of triangles) and remember immediate neighbors for further calculation

        std::vector< std::set< unsigned int > > vertex_neighbors;

        std::map<unsigned int, unsigned int> oddValence; // keep track of the 3rd vertex in case of ordinary odd vertex



        // II) Then, compute the positions for the even vertices: (make sure that you handle the boundaries correctly)-> EVEN MASK

        for (unsigned int v = 0; v < _vertexPositions.size(); v++) {

            //TODO: Compute the coordinates for even vertices - check both the cases - ordinary and extraordinary

            int n = trianglesOnEdge[v].size();
            float alpha_n = (40.0 - pow(3.0 + 2.0 * cos(2.f * M_PI / n), 2)) / 64.0;
            ////Warren
            //if (n == 3) {
            //    beta = 3.0f / 16.0f;
            //}
            //else if (n > 3) {
            //    beta = 3.0f / ((float)n * 8.0f);
            //}

            newVertices[v] *= (1 - alpha_n);

            // If an edge is found only once for a vertice it means this vertice is inside a "open" mesh, i.e. an extraordinary mesh
            bool ordinary = true;
            for (auto it = trianglesOnEdge[v].begin(); it != trianglesOnEdge[v].end(); ++it) {
                // it->first is and Edge object, so it->first.a returns the a vertex of the Edge 
                //and it->first.b returns the b vertex
                unsigned int neighbor_vertex = (it->first.a == v) ? it->first.b : it->first.a;
                if (ordinary) {
                    if (it->second == 2) {
                        newVertices[v] += _vertexPositions[neighbor_vertex] * alpha_n / (float)n;
                    }
                    else {
                        ordinary = false;
                        newVertices[v] = _vertexPositions[v] * 3.f / 4.f + _vertexPositions[neighbor_vertex] / 8.f;
                    }
                }
                else {
                    if (it->second != 2)
                        newVertices[v] += _vertexPositions[neighbor_vertex] / 8.f;
                }
            }
        }


        // III) Then, compute the odd vertices: 
        //se già c'è perche è un edge condiviso ok, altrimenti lo creo il nuovo vertice
        for (unsigned int tIt = 0; tIt < _triangleIndices.size(); ++tIt) {
            unsigned int a = _triangleIndices[tIt][0];
            unsigned int b = _triangleIndices[tIt][1];
            unsigned int c = _triangleIndices[tIt][2];


            //TODO: Update odd vertices ->ODD MASK

            Edge Eab(a, b);
            unsigned int oddVertexOnEdgeEab = 0;
            if (newVertexOnEdge.find(Eab) == newVertexOnEdge.end()) { //we always enter here the first time of the edge 
                //if the vertex doesnt exist, i create it in the middle of the edge to use it in the next iter
                newVertices.push_back((_vertexPositions[a] + _vertexPositions[b]) / 2.f); //it is okay also for the extraordinary case
                oddVertexOnEdgeEab = newVertices.size() - 1; //index of the new vertex
                newVertexOnEdge[Eab] = oddVertexOnEdgeEab;
                oddValence[oddVertexOnEdgeEab] = c;  // keep track of the 3rd vertex if we face an ordinary case

            }
            else { //if it already exixts, update the weight of vertices -> we are in an ordinary case because we meet the edge twice
                oddVertexOnEdgeEab = newVertexOnEdge[Eab];
                newVertices[oddVertexOnEdgeEab] = newVertices[oddVertexOnEdgeEab] * 3.f / 4.f +   // 1/2 * 3/4 = 3/8 (a and b)
                    _vertexPositions[oddValence[oddVertexOnEdgeEab]] / 8.f +  //1/8 
                    _vertexPositions[c] / 8.f; //1/8 
            }



            Edge Ebc(b, c);
            unsigned int oddVertexOnEdgeEbc = 0;
            if (newVertexOnEdge.find(Ebc) == newVertexOnEdge.end()) { //se da true vuol dire che non è contenuto nel vettore
                //se non esiste quell edge, allora creo il vertice al suo centro
                newVertices.push_back((_vertexPositions[b] + _vertexPositions[c]) / 2.f);
                //cerca il nuovo vertice al centro del nuovo edge, da usare nella prossima iterazione
                oddVertexOnEdgeEbc = newVertices.size() - 1; //index of the new vertex
                newVertexOnEdge[Ebc] = oddVertexOnEdgeEbc;
                oddValence[oddVertexOnEdgeEbc] = a;
                // keep track of the 3rd vertex if we face an ordinary case

            }
            else { //if it already exixts, update the weight of vertices
                oddVertexOnEdgeEbc = newVertexOnEdge[Ebc];
                newVertices[oddVertexOnEdgeEbc] = newVertices[oddVertexOnEdgeEbc] * 3.f / 4.f +   // 1/2 * 3/4 = 3/8 (a and b)
                    _vertexPositions[oddValence[oddVertexOnEdgeEbc]] / 8.f +  //1/8 
                    _vertexPositions[a] / 8.f; //1/8 
            }


            Edge Eca(c, a);
            unsigned int oddVertexOnEdgeEca = 0;
            if (newVertexOnEdge.find(Eca) == newVertexOnEdge.end()) {
                //se non esiste quell edge, allora creo il vertice al suo centro
                newVertices.push_back((_vertexPositions[c] + _vertexPositions[a]) / 2.f);
                //cerca il nuovo vertice al centro del nuovo edge, da usare nella prossima iterazione
                oddVertexOnEdgeEca = newVertices.size() - 1; //index of the new vertex
                newVertexOnEdge[Eca] = oddVertexOnEdgeEca;
                oddValence[oddVertexOnEdgeEca] = b;
                // keep track of the 3rd vertex if we face an ordinary case

            }
            else { //if it already exixts, update the weight of vertices
                oddVertexOnEdgeEca = newVertexOnEdge[Eca];
                newVertices[oddVertexOnEdgeEca] = newVertices[oddVertexOnEdgeEca] * 3.f / 4.f +   // 1/2 * 3/4 = 3/8 (a and b)
                    _vertexPositions[oddValence[oddVertexOnEdgeEca]] / 8.f +  //1/8 
                    _vertexPositions[b] / 8.f; //1/8 
            }


            // set new triangles :
            newTriangles.push_back(glm::uvec3(a, oddVertexOnEdgeEab, oddVertexOnEdgeEca));
            newTriangles.push_back(glm::uvec3(oddVertexOnEdgeEab, b, oddVertexOnEdgeEbc));
            newTriangles.push_back(glm::uvec3(oddVertexOnEdgeEca, oddVertexOnEdgeEbc, c));
            newTriangles.push_back(glm::uvec3(oddVertexOnEdgeEab, oddVertexOnEdgeEbc, oddVertexOnEdgeEca));

        }


        // after that:
        _triangleIndices = newTriangles;
        _vertexPositions = newVertices;

        recomputePerVertexNormals();
        recomputePerVertexTextureCoordinates();
    }


    void addNoise() {
        if (_originalVertexPositions.size() == 0) {
            _originalVertexPositions = _vertexPositions;
        } //the first time I copy the orginal vertex pos

        for (unsigned int v = 0; v < _vertexPositions.size(); ++v) {

            _vertexPositions[v] += ((rand() % 100) - 50) * 0.0001;

        }


        recomputePerVertexNormals();
        recomputePerVertexTextureCoordinates();
        _noisyVertexPositions = _vertexPositions;
    }


    float calculateMSE() {
        float mse = 0.f;

        for (unsigned int v = 0; v < _vertexPositions.size(); ++v) {

            glm::vec3 diff = _originalVertexPositions[v] - _vertexPositions[v];
            mse += glm::dot(diff, diff);

        }
        return mse / _vertexPositions.size();
    }


    float findOptimalSmoothingParameters() {
        float bestMSE = std::numeric_limits<float>::max();
        float optimalLambda = 0.f;
        float optimalmu = 0.f;

        // Try different smoothing parameters
        for (float lambda = 0.2f; lambda <= .9f; lambda += 0.1f) {

            for (float mu = lambda + 0.05f; mu <= lambda + 0.25f; mu += 0.05f) {

                // Apply smoothing with the current parameter
                taubinSmoothing(lambda, mu);

                // Calculate MSE between the smoothed mesh and the original mesh
                float mse = calculateMSE();

                // If the MSE is better than the current best, update the best MSE and the optimal smoothing parameter
                if (mse < bestMSE) {

                    bestMSE = mse;
                    optimalLambda = lambda;
                    optimalmu = mu;
                }

                std::cout << "lambda: " << lambda << std::endl;
                std::cout << "mu: " << mu << std::endl;
                std::cout << "MSE: " << mse << std::endl;

                // Reset the mesh to the original state for the next iteration
                _vertexPositions = _noisyVertexPositions;

                recomputePerVertexNormals();
                recomputePerVertexTextureCoordinates();
            }
        }

        std::cout << "OPTIlambda: " << optimalLambda << std::endl;
        std::cout << "OPTImu: " << optimalmu << std::endl;
        std::cout << "OPTIMSE: " << bestMSE << std::endl;

        return optimalLambda, optimalmu;
    }



    void smoothingtaubin2(float lambda, float mu) {


        struct Edge {
            unsigned int a, b;
            Edge(unsigned int c, unsigned int d) : a(std::min<unsigned int>(c, d)), b(std::max<unsigned int>(c, d)) {}
            bool operator < (Edge const& o) const { return a < o.a || (a == o.a && b < o.b); }
            bool operator == (Edge const& o) const { return a == o.a && b == o.b; }
        };

        std::vector< std::map<Edge, unsigned int>> trianglesOnEdge(_vertexPositions.size()); ; // this will be useful to find out if an edge is boundary or not
        std::vector< std::set< unsigned int > > neighboringVertices(_vertexPositions.size()); // this will be used to store the adjacent vertices, i.e., neighboringVertices[i]

        for (unsigned int tIt = 0; tIt < _triangleIndices.size(); ++tIt) {
            unsigned int a = _triangleIndices[tIt][0];
            unsigned int b = _triangleIndices[tIt][1];
            unsigned int c = _triangleIndices[tIt][2];

            //TODO: Remember the faces shared by the edge

            Edge Eab(a, b);
            Edge Ebc(b, c);
            Edge Eac(a, c);

            // a
            trianglesOnEdge[a][Eab] = (trianglesOnEdge[a].find(Eab) == trianglesOnEdge[a].end()) ? 1 : 2;
            trianglesOnEdge[a][Eac] = (trianglesOnEdge[a].find(Eac) == trianglesOnEdge[a].end()) ? 1 : 2;
            // b
            trianglesOnEdge[b][Eab] = (trianglesOnEdge[b].find(Eab) == trianglesOnEdge[b].end()) ? 1 : 2;
            trianglesOnEdge[b][Ebc] = (trianglesOnEdge[b].find(Ebc) == trianglesOnEdge[b].end()) ? 1 : 2;
            // c
            trianglesOnEdge[c][Eac] = (trianglesOnEdge[c].find(Eac) == trianglesOnEdge[c].end()) ? 1 : 2;
            trianglesOnEdge[c][Ebc] = (trianglesOnEdge[c].find(Ebc) == trianglesOnEdge[c].end()) ? 1 : 2;



            //to store the neighbors

            neighboringVertices[a].insert(b); //quindi sul vicino di a sono aggiunti b e c
            neighboringVertices[a].insert(c);
            neighboringVertices[b].insert(a);
            neighboringVertices[b].insert(c);
            neighboringVertices[c].insert(a);
            neighboringVertices[c].insert(b);


        }



        //NON SERVE MEMORIZZARE I VICINI PERCHE 
        //1) CALCOLO vector avarage-> vectorAvg= sum(wij*(vj-vi)
        //wij = 1 / #neighbours;
        //2) CALCOOLO Vi aggiornato-> vNew= vOld+ lambda*vectorAvg
        //(delta tra 0 e 1)
        //
        //RIPETO LA STESSA COSA MA IN 2 CAMBIO lambda CON mu che è negativo e maggiore di lambda in abs
        int N = 25;
        for (int n = 0; n < N; ++n) {

            for (unsigned int v = 0; v < _vertexPositions.size(); ++v) {

                float wij = 0.f;
                float lij = 0.f;

                glm::vec3 vectorAvg(0.0f);

                for (auto it = neighboringVertices[v].begin(); it != neighboringVertices[v].end(); ++it) {

                    lij += neighboringVertices[*it].size();


                    vectorAvg += (_vertexPositions[*it] - _vertexPositions[v]);

                }
                vectorAvg = (1.f / lij) * vectorAvg;

                _vertexPositions[v] += lambda * vectorAvg;

            }

            recomputePerVertexNormals();
            recomputePerVertexTextureCoordinates();



            for (unsigned int v = 0; v < _vertexPositions.size(); ++v) {

                float wij = 0.f;
                float lij = 0.f;

                glm::vec3 vectorAvg(0.0f);

                for (auto it = neighboringVertices[v].begin(); it != neighboringVertices[v].end(); ++it) {

                    lij += neighboringVertices[*it].size();


                    vectorAvg += (_vertexPositions[*it] - _vertexPositions[v]);

                }
                vectorAvg = (1.f / lij) * vectorAvg;

                _vertexPositions[v] += -mu * vectorAvg;


            }


            recomputePerVertexNormals();
            recomputePerVertexTextureCoordinates();
        }
    }




    void taubinSmoothing(float lambda, float mu) {

        struct Edge {
            unsigned int a, b;
            Edge(unsigned int c, unsigned int d) : a(std::min<unsigned int>(c, d)), b(std::max<unsigned int>(c, d)) {}
            bool operator < (Edge const& o) const { return a < o.a || (a == o.a && b < o.b); }
            bool operator == (Edge const& o) const { return a == o.a && b == o.b; }
        };
        std::vector< std::map<Edge, unsigned int>> trianglesOnEdge(_vertexPositions.size()); ; // this will be useful to find out if an edge is boundary or not
        std::vector< std::set< unsigned int > > neighboringVertices(_vertexPositions.size()); // this will be used to store the adjacent vertices, i.e., neighboringVertices[i]

        for (unsigned int tIt = 0; tIt < _triangleIndices.size(); ++tIt) {
            unsigned int a = _triangleIndices[tIt][0];
            unsigned int b = _triangleIndices[tIt][1];
            unsigned int c = _triangleIndices[tIt][2];

            //TODO: Remember the faces shared by the edge

            Edge Eab(a, b);
            Edge Ebc(b, c);
            Edge Eac(a, c);

            // a
            trianglesOnEdge[a][Eab] = (trianglesOnEdge[a].find(Eab) == trianglesOnEdge[a].end()) ? 1 : 2;
            trianglesOnEdge[a][Eac] = (trianglesOnEdge[a].find(Eac) == trianglesOnEdge[a].end()) ? 1 : 2;
            // b
            trianglesOnEdge[b][Eab] = (trianglesOnEdge[b].find(Eab) == trianglesOnEdge[b].end()) ? 1 : 2;
            trianglesOnEdge[b][Ebc] = (trianglesOnEdge[b].find(Ebc) == trianglesOnEdge[b].end()) ? 1 : 2;
            // c
            trianglesOnEdge[c][Eac] = (trianglesOnEdge[c].find(Eac) == trianglesOnEdge[c].end()) ? 1 : 2;
            trianglesOnEdge[c][Ebc] = (trianglesOnEdge[c].find(Ebc) == trianglesOnEdge[c].end()) ? 1 : 2;



            //to store the neighbors

            neighboringVertices[a].insert(b); //quindi sul vicino di a sono aggiunti b e c
            neighboringVertices[a].insert(c);
            neighboringVertices[b].insert(a);
            neighboringVertices[b].insert(c);
            neighboringVertices[c].insert(a);
            neighboringVertices[c].insert(b);


        }

        int N = 10;
        for (int n = 0; n < N; ++n) {

            for (unsigned int v = 0; v < _vertexPositions.size(); ++v) {
                for (auto it = trianglesOnEdge[v].begin(); it != trianglesOnEdge[v].end(); ++it) {
                    // it->first is and Edge object
                    // it->second the number of triangles
                  
                    if (it->second == 2) {

                        glm::vec3 delta(0.0f);
                        float wij = 0.f;

                        for (auto i = neighboringVertices[v].begin(); i != neighboringVertices[v].end(); ++i) {

                            wij = 1.f / neighboringVertices[*i].size();
                            delta += wij * (_vertexPositions[*i] - _vertexPositions[v]);

                        }
                        // shrink step
                        _vertexPositions[v] += (lambda * delta);
                    }

                }
            }

           recomputePerVertexNormals();
           recomputePerVertexTextureCoordinates();

            _lambdaVertexNormals = _vertexNormals;


            for (unsigned int v = 0; v < _vertexPositions.size(); ++v) {
              
                        glm::vec3 delta(0.0f);
                        float wij = 0.f;

                        for (auto i = neighboringVertices[v].begin(); i != neighboringVertices[v].end(); ++i) {

                            wij = 1.f / neighboringVertices[*i].size();
                            delta += wij * (_vertexPositions[*i] - _vertexPositions[v]);
                  
                        }
                        // shrink step
                        _vertexPositions[v] += (-mu * delta);
               
            }

            recomputePerVertexNormals();
            recomputePerVertexTextureCoordinates();




        }



        int sum = 0;
        int notIn = 0;

        for (unsigned int v = 0; v < _vertexPositions.size(); ++v) {

            float dotProduct = glm::dot(_vertexNormals[v], _lambdaVertexNormals[v]);

            // confronta il prodotto interno con 0.1
            if (dotProduct < 0.001) {
                sum += 1;
            }
            else {
                notIn += 1;
            }

        }


    }



    void laplacianSmoothing(float lambda) {

        struct Edge {
            unsigned int a, b;
            Edge(unsigned int c, unsigned int d) : a(std::min<unsigned int>(c, d)), b(std::max<unsigned int>(c, d)) {}
            bool operator < (Edge const& o) const { return a < o.a || (a == o.a && b < o.b); }
            bool operator == (Edge const& o) const { return a == o.a && b == o.b; }
        };
        std::vector< std::map<Edge, unsigned int>> trianglesOnEdge(_vertexPositions.size()); ; // this will be useful to find out if an edge is boundary or not
        std::vector< std::set< unsigned int > > neighboringVertices(_vertexPositions.size()); // this will be used to store the adjacent vertices, i.e., neighboringVertices[i]

        for (unsigned int tIt = 0; tIt < _triangleIndices.size(); ++tIt) {
            unsigned int a = _triangleIndices[tIt][0];
            unsigned int b = _triangleIndices[tIt][1];
            unsigned int c = _triangleIndices[tIt][2];

            //TODO: Remember the faces shared by the edge

            Edge Eab(a, b);
            Edge Ebc(b, c);
            Edge Eac(a, c);

            // a
            trianglesOnEdge[a][Eab] = (trianglesOnEdge[a].find(Eab) == trianglesOnEdge[a].end()) ? 1 : 2;
            trianglesOnEdge[a][Eac] = (trianglesOnEdge[a].find(Eac) == trianglesOnEdge[a].end()) ? 1 : 2;
            // b
            trianglesOnEdge[b][Eab] = (trianglesOnEdge[b].find(Eab) == trianglesOnEdge[b].end()) ? 1 : 2;
            trianglesOnEdge[b][Ebc] = (trianglesOnEdge[b].find(Ebc) == trianglesOnEdge[b].end()) ? 1 : 2;
            // c
            trianglesOnEdge[c][Eac] = (trianglesOnEdge[c].find(Eac) == trianglesOnEdge[c].end()) ? 1 : 2;
            trianglesOnEdge[c][Ebc] = (trianglesOnEdge[c].find(Ebc) == trianglesOnEdge[c].end()) ? 1 : 2;



            //to store the neighbors

            neighboringVertices[a].insert(b); //quindi sul vicino di a sono aggiunti b e c
            neighboringVertices[a].insert(c);
            neighboringVertices[b].insert(a);
            neighboringVertices[b].insert(c);
            neighboringVertices[c].insert(a);
            neighboringVertices[c].insert(b);


        }


        float Ni = 0.f;
        glm::vec3 neigh(0.f);
        glm::vec3 L(0.f);


        int N = 25;
        for (int n = 0; n < N; ++n) {

            for (unsigned int v = 0; v < _vertexPositions.size(); ++v) {


                for (auto it = trianglesOnEdge[v].begin(); it != trianglesOnEdge[v].end(); ++it) {

                    if (it->second == 2) {
                        for (auto i = neighboringVertices[v].begin(); i != neighboringVertices[v].end(); ++i) {
                            L += _vertexPositions[*i] - _vertexPositions[v];
                        }
                        Ni = neighboringVertices[v].size();
                        L /= Ni;

                        _vertexPositions[v] += lambda * L;
                    }
                }
            }

            recomputePerVertexNormals();
            recomputePerVertexTextureCoordinates();
        }
    }



    void gaussianSmoothing(float sigma) {

        struct Edge {
            unsigned int a, b;
            Edge(unsigned int c, unsigned int d) : a(std::min<unsigned int>(c, d)), b(std::max<unsigned int>(c, d)) {}
            bool operator < (Edge const& o) const { return a < o.a || (a == o.a && b < o.b); }
            bool operator == (Edge const& o) const { return a == o.a && b == o.b; }
        };
        std::vector< std::map<Edge, unsigned int>> trianglesOnEdge(_vertexPositions.size()); ; // this will be useful to find out if an edge is boundary or not
        std::vector< std::set< unsigned int > > neighboringVertices(_vertexPositions.size()); // this will be used to store the adjacent vertices, i.e., neighboringVertices[i]

        for (unsigned int tIt = 0; tIt < _triangleIndices.size(); ++tIt) {
            unsigned int a = _triangleIndices[tIt][0];
            unsigned int b = _triangleIndices[tIt][1];
            unsigned int c = _triangleIndices[tIt][2];

            //TODO: Remember the faces shared by the edge

            Edge Eab(a, b);
            Edge Ebc(b, c);
            Edge Eac(a, c);

            // a
            trianglesOnEdge[a][Eab] = (trianglesOnEdge[a].find(Eab) == trianglesOnEdge[a].end()) ? 1 : 2;
            trianglesOnEdge[a][Eac] = (trianglesOnEdge[a].find(Eac) == trianglesOnEdge[a].end()) ? 1 : 2;
            // b
            trianglesOnEdge[b][Eab] = (trianglesOnEdge[b].find(Eab) == trianglesOnEdge[b].end()) ? 1 : 2;
            trianglesOnEdge[b][Ebc] = (trianglesOnEdge[b].find(Ebc) == trianglesOnEdge[b].end()) ? 1 : 2;
            // c
            trianglesOnEdge[c][Eac] = (trianglesOnEdge[c].find(Eac) == trianglesOnEdge[c].end()) ? 1 : 2;
            trianglesOnEdge[c][Ebc] = (trianglesOnEdge[c].find(Ebc) == trianglesOnEdge[c].end()) ? 1 : 2;



            //to store the neighbors

            neighboringVertices[a].insert(b); //quindi sul vicino di a sono aggiunti b e c
            neighboringVertices[a].insert(c);
            neighboringVertices[b].insert(a);
            neighboringVertices[b].insert(c);
            neighboringVertices[c].insert(a);
            neighboringVertices[c].insert(b);


        }


        
        glm::vec3 L(0.f);


        int N = 25;
        for (int n = 0; n < N; ++n) {

            for (unsigned int v = 0; v < _vertexPositions.size(); ++v) {


                for (auto it = trianglesOnEdge[v].begin(); it != trianglesOnEdge[v].end(); ++it) {
                   
                    if (it->second == 2) {

                        float sum = 0.f;
                        glm::vec3 newPosition(0.0f);

                        for (auto i = neighboringVertices[v].begin(); i != neighboringVertices[v].end(); ++i) {
                            L = _vertexPositions[v] - _vertexPositions[*i];
                            float wij = exp(-(glm::dot(L, L)) / (2.0f * sigma * sigma));

                            newPosition += wij * _vertexPositions[*i];
                            sum += wij;

                        }
                         newPosition /= sum;
                         _vertexPositions[v] = newPosition;
                    }
                }
            }

            recomputePerVertexNormals();
            recomputePerVertexTextureCoordinates();
        }
    }



    void cotanLaplacianSmoothing(float lambda) {
        // Declare new vertices and new triangles. Initialize the new positions for the even vertices with (0,0,0):

        struct Edge {
            unsigned int a, b;
            Edge(unsigned int c, unsigned int d) : a(std::min<unsigned int>(c, d)), b(std::max<unsigned int>(c, d)) {}
            bool operator < (Edge const& o) const { return a < o.a || (a == o.a && b < o.b); }
            bool operator == (Edge const& o) const { return a == o.a && b == o.b; }
        };
        std::vector< std::map<Edge, unsigned int>> trianglesOnEdge(_vertexPositions.size()); ; // this will be useful to find out if an edge is boundary or not
        std::vector< std::set< unsigned int > > neighboringVertices(_vertexPositions.size()); // this will be used to store the adjacent vertices, i.e., neighboringVertices[i]

        for (unsigned int tIt = 0; tIt < _triangleIndices.size(); ++tIt) {
            unsigned int a = _triangleIndices[tIt][0];
            unsigned int b = _triangleIndices[tIt][1];
            unsigned int c = _triangleIndices[tIt][2];

            //TODO: Remember the faces shared by the edge

            Edge Eab(a, b);
            Edge Ebc(b, c);
            Edge Eac(a, c);

            // a
            trianglesOnEdge[a][Eab] = (trianglesOnEdge[a].find(Eab) == trianglesOnEdge[a].end()) ? 1 : 2;
            trianglesOnEdge[a][Eac] = (trianglesOnEdge[a].find(Eac) == trianglesOnEdge[a].end()) ? 1 : 2;
            // b
            trianglesOnEdge[b][Eab] = (trianglesOnEdge[b].find(Eab) == trianglesOnEdge[b].end()) ? 1 : 2;
            trianglesOnEdge[b][Ebc] = (trianglesOnEdge[b].find(Ebc) == trianglesOnEdge[b].end()) ? 1 : 2;
            // c
            trianglesOnEdge[c][Eac] = (trianglesOnEdge[c].find(Eac) == trianglesOnEdge[c].end()) ? 1 : 2;
            trianglesOnEdge[c][Ebc] = (trianglesOnEdge[c].find(Ebc) == trianglesOnEdge[c].end()) ? 1 : 2;



            //to store the neighbors

            neighboringVertices[a].insert(b); //quindi sul vicino di a sono aggiunti b e c
            neighboringVertices[a].insert(c);
            neighboringVertices[b].insert(a);
            neighboringVertices[b].insert(c);
            neighboringVertices[c].insert(a);
            neighboringVertices[c].insert(b);



        }


        float Ni = 0.f;

        int N = 5;
        for (int n = 0; n < N; ++n) {

            for (unsigned int v = 0; v < _vertexPositions.size(); ++v) {
                glm::vec3 L(0.f);
                glm::vec3 delta(0.f);
                float sum = 0.f;
                float area = 0.f;

                for (auto it = trianglesOnEdge[v].begin(); it != trianglesOnEdge[v].end(); ++it) {

                    if (it->second == 2) {
                        for (auto i = neighboringVertices[v].begin(); i != neighboringVertices[v].end(); ++i) {


                            unsigned int i_next = (std::next(i) == neighboringVertices[v].end()) ? *(neighboringVertices[v].begin()) : *(std::next(i));
                            unsigned int i_prev = (i == neighboringVertices[v].begin()) ? *(std::prev(neighboringVertices[v].end())) : *(std::prev(i));

                            glm::vec3 vi = _vertexPositions[v];
                            glm::vec3 vj = _vertexPositions[*i]; //i,j
                            glm::vec3 viprev = _vertexPositions[i_prev]; //i,j-1
                            glm::vec3 vinext = _vertexPositions[i_next]; //i,j

                            float e1 = glm::length(vi - vj);
                            float e2 = glm::length(vi - viprev);
                            float e3 = glm::length(viprev - vj);
                           
                            area = 0.5f * glm::length(glm::cross(vj - vi, viprev - vi));

                            // NOTE: cos(alpha) = (a^2.b^2  - c^2) / (2.a.b)
                           // with a, b, c the lengths of of the sides of the triangle and (a, b)
                           // forming the angle alpha.

                            float cos_alpha = fabs((e3 * e3 + e2 * e2 - e1 * e1) / (2.0f * e3 * e2));

                            float e4 = glm::length(vi - vinext);
                            float e5 = glm::length(vinext - vj);
                            float cos_beta = fabs((e4 * e4 + e5 * e5 - e1 * e1) / (2.0f * e4 * e5));
                            // cot(x) = cos(x)/sin(x)
                            // cos(x)^2 + sin(x)^2 = 1
                            // then sin(x) = sqrt(1-cos(x)^2)

                            float cotan1 = cos_alpha / sqrt(1.0f - cos_alpha * cos_alpha);
                            float cotan2 = cos_beta / sqrt(1.0f - cos_beta * cos_beta);
                            float wij = (cotan1 + cotan2) * 0.5f;

                            if (isnan(wij)) {
                                wij = 0.0f;
                            }
                            // cotangent value close to 0.0f.
                            // as cotan approaches infinity close to zero we clamp
                            // higher values.333333333333333<
                            const float eps = 1e-6f;
                            const float cotan_max = cos(eps) / sin(eps);
                            if (wij >= cotan_max) {
                                wij = cotan_max;
                            }

                            delta += wij * (_vertexPositions[*i]);
                            sum += wij;

                        }
                        
                        float lambda = 0.3f;
                        _vertexPositions[v] = (lambda * (delta/sum)) + (1.f - lambda)* _vertexPositions[v];

                    } 

                }
            }


            recomputePerVertexNormals();
            recomputePerVertexTextureCoordinates();
        }

    }





    void subdivide() {

        //subdivideLinear();
        //for (int i = 0; i < 2; i++)   
            subdivideLoop();

    }

private:
    std::vector<glm::vec3> _vertexPositions;
    std::vector<glm::vec3> _originalVertexPositions;
    std::vector<glm::vec3> _noisyVertexPositions;
    std::vector<glm::vec3> _vertexNormals;
    std::vector<glm::vec3> _lambdaVertexNormals;
    std::vector<glm::vec2> _vertexTexCoords;
    std::vector<glm::uvec3> _triangleIndices;

    GLuint _vao = 0;
    GLuint _posVbo = 0;
    GLuint _normalVbo = 0;
    GLuint _texCoordVbo = 0;
    GLuint _ibo = 0;
};

// utility: loader
void loadOFF(const std::string& filename, std::shared_ptr<Mesh> meshPtr);

#endif  // MESH_H
