///////////////////////////////////////////////////
//
//  Hamish Carr
//  September, 2020
//
//  ------------------------
//  DirectedEdgeSurface.cpp
//  ------------------------
//  
//  Base code for rendering assignments.
//
//  Minimalist (non-optimised) code for reading and 
//  rendering an object file
//  
//  We will make some hard assumptions about input file
//  quality. We will not check for manifoldness or 
//  normal direction, &c.  And if it doesn't work on 
//  all object files, that's fine.
//
//  While I could set it up to use QImage for textures,
//  I want this code to be reusable without Qt, so I 
//  shall make a hard assumption that textures are in 
//  ASCII PPM and use my own code to read them
//  
///////////////////////////////////////////////////

// include the header file
#include "DirectedEdgeSurface.h"

// include the C++ standard libraries we want
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

// include the Cartesian 3- vector class
#include "Cartesian3.h"
#include "SphereVertices.h"

#define MAXIMUM_LINE_LENGTH 1024

// constructor will initialise to safe values
DirectedEdgeSurface::DirectedEdgeSurface()
    : centreOfGravity(0.0,0.0,0.0)
    { // DirectedEdgeSurface()
    // force arrays to size 0
    vertices.resize(0);
    normals.resize(0);
	firstDirectedEdge.resize(0);
	faceVertices.resize(0);
	otherHalf.resize(0);
    } // DirectedEdgeSurface()

// read routine returns true on success, failure otherwise
bool DirectedEdgeSurface::ReadObjectStream(std::istream &geometryStream)
    { // ReadObjectStream()
    
	// create a read buffer
	char readBuffer[MAXIMUM_LINE_LENGTH];

	while (true)
	{ // not eof
	// token for identifying meaning of line
		std::string token;

		// character to read
		geometryStream >> token;

		// check for eof() in case we've run out
		if (geometryStream.eof())
			break;

		// otherwise, switch on the token we read
		if (token == "#")
		{ // comment 
		// read and discard the line
			geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
		} // comment
		else if (token == "Vertex")
		{ // vertex
		// variables for the read
			unsigned int vertexID;
			geometryStream >> vertexID;
			// it has to be next valid 0-based ID, so
			// reject line if it isn't
			if (vertexID != vertices.size())
			{ // bad vertex ID
			// read and discard the line
				geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
			} // bad vertex ID				

		// read in the new vertex position
			Cartesian3 newVertex;
			geometryStream >> newVertex;
			HE_Vertex new_vert;
			new_vert.pos = newVertex;
			// and add it to the vertices
			vertices.push_back(new_vert);
		} // vertex
		else if (token == "Normal")
		{ // normal
		// variables for the read
			unsigned int normalID;
			geometryStream >> normalID;
			// it has to be next valid 0-based ID, so
			// reject line if it isn't
			if (normalID != normals.size())
			{ // bad ID
			// read and discard the line
				geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
			} // bad ID				

		// read in the new normal
			Cartesian3 newNormal;
			geometryStream >> newNormal;

			// and add it to the vertices
			normals.push_back(newNormal);
		} // normal
		else if (token == "FirstDirectedEdge")
		{ // first directed edge
		// variables for the read
			unsigned int FDEID;
			geometryStream >> FDEID;
			// it has to be next valid 0-based ID, so
			// reject line if it isn't
			if (vertices[FDEID].FDedge != -1)
			{ // bad ID
			// read and discard the line
				geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
			} // bad ID				

		// read in the new FDE
			unsigned int newFDE;
			geometryStream >> newFDE;

			// and add it to the vertices
			vertices[FDEID].FDedge = newFDE;
		} // first directed edge
		else if (token == "Face")
		{ // face
		// variables for the read
			unsigned int faceID;
			geometryStream >> faceID;
			// it has to be next valid 0-based ID, so
			// reject line if it isn't
			if (faceID != faceVertices.size() / 3)
			{ // bad face ID
			// read and discard the line
				geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
			} // bad face ID				

		// read in the new face vertex (3 times)
			unsigned int newFaceVertex;
			geometryStream >> newFaceVertex;
			faceVertices.push_back(newFaceVertex);
			geometryStream >> newFaceVertex;
			faceVertices.push_back(newFaceVertex);
			geometryStream >> newFaceVertex;
			faceVertices.push_back(newFaceVertex);
		} // face
		else if (token == "OtherHalf")
		{ // other half
		// variables for the read
			unsigned int otherHalfID;
			geometryStream >> otherHalfID;
			// it has to be next valid 0-based ID, so
			// reject line if it isn't
			if (otherHalfID != halfEdges.size())
			{ // bad ID
			// read and discard the line
				geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
			} // bad ID				

		// read in the new face vertex (3 times)
			unsigned int newOtherHalf;
			geometryStream >> newOtherHalf;
			HE_halfEdge newOther;
			newOther.otherHalf = newOtherHalf;
			newOther.face = otherHalfID / 3;
			newOther.begin = faceVertices[otherHalfID];
			//judge next half-edge index ,order: 0-1 1-2 2-0
			if (otherHalfID % 3 == 0)
			{
				newOther.next = otherHalfID + 1;
				newOther.prev = otherHalfID + 2;
			}
			else if (otherHalfID % 3 == 1)
			{
				newOther.next = otherHalfID + 1;
				newOther.prev = otherHalfID - 1;
			}
			else if (otherHalfID % 3 == 2)
			{
				newOther.next = otherHalfID - 2;
				newOther.prev = otherHalfID - 1;
			}
			halfEdges.push_back(newOther);

		} // other half
	} // not eof

// compute centre of gravity
// note that very large files may have numerical problems with this
	centreOfGravity = Cartesian3(0.0, 0.0, 0.0);

	// if there are any vertices at all
	if (vertices.size() != 0)
	{ // non-empty vertex set
	// sum up all of the vertex positions
		for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
			centreOfGravity = centreOfGravity + vertices[vertex].pos;

		// and divide through by the number to get the average position
		// also known as the barycentre
		centreOfGravity = centreOfGravity / vertices.size();

		// start with 0 radius
		objectSize = 0.0;

		// now compute the largest distance from the origin to a vertex
		for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
		{ // per vertex
		// compute the distance from the barycentre
			float distance = (vertices[vertex].pos - centreOfGravity).length();

			// now test for maximality
			if (distance > objectSize)
				objectSize = distance;
		} // per vertex
	} // non-empty vertex set

// return a success code
	return true;
    } // ReadObjectStream()


// write routine
void DirectedEdgeSurface::WriteObjectStream(std::ostream *geometryStream)
    { // WriteObjectStream()
	*geometryStream << "#" << std::endl;
	*geometryStream << "# Loop Subdivision" << std::endl;
	*geometryStream << "# Chenxin Jiang" << std::endl;
	*geometryStream << "# 201614769" << std::endl;

	// output the vertices
	for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
		*geometryStream << "Vertex " << vertex << " " << std::fixed << vertices[vertex].pos << std::endl;

	// and the normal vectors
	for (unsigned int normal = 0; normal < normals.size(); normal++)
		*geometryStream << "Normal " << normal << " " << std::fixed << normals[normal] << std::endl;

	// and the first directed edges
	for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
		*geometryStream << "FirstDirectedEdge " << vertex << " " << std::fixed << vertices[vertex].FDedge << std::endl;

	// and the faces - increment is taken care of internally
	for (unsigned int face = 0; face < faceVertices.size(); )
	{ // per face
		*geometryStream << "Face " << face / 3 << " ";

		// read in three vertices
		*geometryStream << faceVertices[face++] << " ";
		*geometryStream << faceVertices[face++] << " ";
		*geometryStream << faceVertices[face++];

		*geometryStream << std::endl;
	} // per face
	*geometryStream << "edge " << "begin " << "face " << "next " << "otherhalf " << "prev " << std::endl;
	for (unsigned int i = 0; i < halfEdges.size(); i++)
		*geometryStream << "edge    " << halfEdges[i].begin << "   " << halfEdges[i].face << "   " << halfEdges[i].next << "     " << halfEdges[i].otherHalf << "      " << halfEdges[i].prev << std::endl;
}// WriteObjectStream()

// routine to render
void DirectedEdgeSurface::Render(RenderParameters *renderParameters)
    { // Render()
    // Ideally, we would apply a global transformation to the object, but sadly that breaks down
    // when we want to scale things, as unless we normalise the normal vectors, we end up affecting
    // the illumination.  Known solutions include:
    // 1.   Normalising the normal vectors
    // 2.   Explicitly dividing the normal vectors by the scale to balance
    // 3.   Scaling only the vertex position (slower, but safer)
    // 4.   Not allowing spatial zoom (note: sniper scopes are a modified projection matrix)
    //
    // Inside a game engine, zoom usually doesn't apply. Normalisation of normal vectors is expensive,
    // so we will choose option 2.  

    // Scale defaults to the zoom setting
    float scale = renderParameters->zoomScale;
    scale /= objectSize;
        
    //  now scale everything
    glScalef(scale, scale, scale);

    // apply the translation to the centre of the object if requested
    glTranslatef(-centreOfGravity.x, -centreOfGravity.y, -centreOfGravity.z);

    // start rendering
    glBegin(GL_TRIANGLES);

	// set colour for pick render - ignored for regular render
	glColor3f(1.0, 1.0, 1.0);

    // loop through the faces
	for (unsigned int face = 0; face < faceVertices.size(); face +=3)
		{ // per face
		// if we want flat normals, compute them here
		if (renderParameters->useFlatNormals)
			{ // flat normals
			// find two vectors along edges of the triangle
			Cartesian3 pq = vertices[faceVertices[face+1]].pos - vertices[faceVertices[face]].pos;
			Cartesian3 pr = vertices[faceVertices[face+2]].pos - vertices[faceVertices[face]].pos;

			// take their cross product and normalise
			Cartesian3 faceNormal = pq.cross(pr).unit();

			// and use it to set the glNormal
			glNormal3f(faceNormal.x * scale, faceNormal.y * scale, faceNormal.z * scale);
			} // flat normals

		// we have made a HARD assumption that we have enough normals
		for (unsigned int vertex = face; vertex < face+3; vertex++)
			{ // per vertex
		
			// if we are using smooth normals
			if (!renderParameters->useFlatNormals)
				// set the normal vector
				glNormal3f
					(
					normals[faceVertices[vertex]].x* scale,
					normals[faceVertices[vertex]].y * scale,
					normals[faceVertices[vertex]].z * scale
					);
			
			// and set the vertex position
			glVertex3f
				(
				vertices[faceVertices[vertex]].pos.x,
				vertices[faceVertices[vertex]].pos.y,
				vertices[faceVertices[vertex]].pos.z
				);

			} // per vertex

		} // per face

    // close off the triangles
    glEnd();
    
    // now we add a second loop to render the vertices if desired
    if (!renderParameters->showVertices)
    	return;

	glDisable(GL_LIGHTING);

	// loop through the vertices
	for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
		{ // per vertex
		// use modelview matrix (not most efficient solution, but quickest to code)
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		glTranslatef(vertices[vertex].pos.x, vertices[vertex].pos.y, vertices[vertex].pos.z);
		glScalef(0.1 * renderParameters->vertexSize, 0.1 * renderParameters->vertexSize, 0.1 * renderParameters->vertexSize);
		renderTriangulatedSphere();
		glPopMatrix();
		} // per vertex 
    
    } // Render()



void DirectedEdgeSurface::updateHalfEdges()
{
	for (int i = 0; i < faceVertices.size() / 3; i++)
	{

		// begin;
		// otherHalf;
		// next;
		// prev;
		// face;

		HE_halfEdge e0, e1, e2;
		e0.begin = faceVertices[i * 3];
		e1.begin = faceVertices[i * 3 + 1];
		e2.begin = faceVertices[i * 3 + 2];

		e0.face = i;
		e1.face = i;
		e2.face = i;

		e0.next = i * 3 + 1;
		e1.next = i * 3 + 2;
		e2.next = i * 3;

		e0.prev = i * 3 + 2;
		e1.prev = i * 3;
		e2.prev = i * 3 + 1;

		halfEdges.push_back(e0);
		halfEdges.push_back(e1);
		halfEdges.push_back(e2);
	}
}

void DirectedEdgeSurface::LoopSubdivision()
{
	//1.create new vertices while upadating position
	for (auto& edge : halfEdges) {
		if (!edge.hasNewVertex) {
			//generate new vertex
			HE_Vertex newVertex;
			newVertex.is_new = true;
			edge.newVertexID = vertices.size();
			halfEdges[edge.otherHalf].newVertexID = vertices.size();
			edge.hasNewVertex = true;
			halfEdges[edge.otherHalf].hasNewVertex = true;


			///Use specific formula to upadate the new vertex position
			newVertex.pos = 0.375 * (vertices[edge.begin].pos + vertices[halfEdges[edge.otherHalf].begin].pos)
				+ 0.125 * (vertices[halfEdges[edge.next].begin].pos + vertices[halfEdges[halfEdges[edge.otherHalf].next].begin].pos);

			vertices.push_back(newVertex);
		}
	}

	//2.update the position of old vertex
	for (int i = 0; i < vertices.size(); i++) {
		if (vertices[i].is_new == false) {
			adjustVertex(i);
		}
		else {
			vertices[i].is_new = false;
		}
	}

	//3.update and add new faces
	//		v2
	//      /\
	//   v6/__\v5
	//    /\v4/\
	// v3/__\/__\ v1
	for (int i = 0; i < halfEdges.size() / 3; i++) {
		//six vertex index
		unsigned int v1 = faceVertices[3 * i];
		unsigned int v2 = faceVertices[3 * i + 1];
		unsigned int v3 = faceVertices[3 * i + 2];
		unsigned int v4 = halfEdges[3 * i].newVertexID;
		unsigned int v5 = halfEdges[3 * i + 1].newVertexID;
		unsigned int v6 = halfEdges[3 * i + 2].newVertexID;

		//edit original face
		faceVertices[3 * i] = v4;
		faceVertices[3 * i + 1] = v5;
		faceVertices[3 * i + 2] = v6;

		//update new faces
		//v1 v5 v4
		faceVertices.push_back(v1);
		faceVertices.push_back(v5);
		faceVertices.push_back(v4);

		//v2 v6 v5
		faceVertices.push_back(v2);
		faceVertices.push_back(v6);
		faceVertices.push_back(v5);

		//v3 v4 v6
		faceVertices.push_back(v3);
		faceVertices.push_back(v4);
		faceVertices.push_back(v6);

	}

	//4.1 rebuild half-edge 
    halfEdges.clear();
	updateHalfEdges();

	//4.2 update other-half
	for (unsigned int i = 0; i < halfEdges.size(); i++) {
		unsigned int start1 = halfEdges[halfEdges[i].prev].begin;
		unsigned int end1 = halfEdges[i].begin;
		for (unsigned int j = 0; j < halfEdges.size(); j++) {
			unsigned int start2 = halfEdges[halfEdges[j].prev].begin;
			unsigned int end2 = halfEdges[j].begin;
			if (start1 == end2 && start2 == end1) {
				halfEdges[i].otherHalf = j;
				break;
			}
		}
	}

	//5.update first directed edge
	for (unsigned int i = 0; i < vertices.size(); i++) {
		for (unsigned int j = 0; j < halfEdges.size(); j++)
		{
			unsigned int start = halfEdges[halfEdges[j].prev].begin;
			if (start == i)
			{
				vertices[i].FDedge = j;
				break;
			}
		}
	}
}

void DirectedEdgeSurface::adjustVertex(int i)
{
    float u;
    Cartesian3 sum = {0, 0, 0};  // Assuming Cartesian3 has 3 float components
    std::vector<HE_Vertex> list1 = AdjacentVertex(i);
    unsigned int degree = list1.size();

    if (degree == 3) {
        u = 0.1875;
    } else {
        u = 3.f / (8.f * degree);  // You might want to adjust this formula
    }

    float precomputed = float(degree) * u;
    for (const auto& vertex : list1) {
        sum = sum + vertex.pos;
    }

    vertices[i].pos = (1.0f - precomputed) * vertices[i].pos + u * sum;
}


std::vector<HE_Vertex> DirectedEdgeSurface::AdjacentVertex(int i)
{
	//get all vector of adjacent vertex
	std::vector<HE_Vertex> adjacent;
	unsigned int FDedge = vertices[i].FDedge;
	unsigned int FDedgeID = FDedge;
	do
	{
		adjacent.push_back(vertices[halfEdges[FDedgeID].begin]);
		FDedgeID = halfEdges[FDedgeID].otherHalf;
		if (halfEdges[FDedgeID].begin == i)
		{
			FDedgeID = halfEdges[FDedgeID].next;
		}
		else { break; }
	} while (FDedgeID != FDedge);
	return adjacent;
}






