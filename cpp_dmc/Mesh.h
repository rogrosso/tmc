#pragma once

#include <iostream>
#include <vector>
#include <array>
#include <map>
#include <set>
#include <algorithm>
#include <iterator>

// project
#include "Vector.h"

namespace homotopy {
	class Mesh {
	public:
		typedef typename Vector Vertex;
		typedef typename Vector Normal;
		typedef typename Vector Point;
		typedef typename std::array<int, 3> Triangle;
		typedef typename std::array<int, 4> Quadrilateral;
		typedef typename std::array<int, 2> Line;
		typedef typename std::set<int> OneRing;
		typedef typename std::array<double, 4> Attribute;
		typedef typename std::array<Vertex, 8> BBox;
	public:
		void resizeVertices(const int size_)
		{
			m_vertices.resize(size_);
			m_normals.resize(size_);
			m_attributes.resize(size_);
		}
		void addVertex(const int pos, Vertex& v, Normal& n)
		{
			m_vertices[pos] = v;
			m_normals[pos] = n;
			/*const double r = 250.0 / 255.0;
			const double g = 235.0 / 255.0;
			const double b = 215.0 / 255.0;
			const double a = 1.0;*/
            const double r = 175.0 / 255.0;
            const double g = 170.0 / 255.0;
            const double b = 170.0 / 255.0;
            const double a = 1.0;
			/*
            const double r = 250.0 / 255.0;
            const double g = 250.0 / 255.0;
            const double b = 250.0 / 255.0;
            const double a = 1.0;
			*/
			m_attributes[pos] = std::array<double,4>{ r,g,b,a };
		}
		void addVertex(Vertex& v, Normal& n)
		{
			m_vertices.push_back(v);
			m_normals.push_back(n);
			//const double r = 250.0 / 255.0;
			//const double g = 235.0 / 255.0;
			//const double b = 215.0 / 255.0;
			//const double a = 1.0;
            /*const double r = 235.0 / 255.0;
            const double g = 235.0 / 255.0;
            const double b = 235.0 / 255.0;
            const double a = 1.0;*/
            const double r = 175.0 / 255.0;
            const double g = 170.0 / 255.0;
            const double b = 170.0 / 255.0;
            const double a = 1.0;
			/*
            const double r = 250.0 / 255.0;
            const double g = 250.0 / 255.0;
            const double b = 250.0 / 255.0;
            const double a = 1.0;
			*/
			m_attributes.push_back({ r,g,b,a });
		}
		void addVertex(Vertex& v, Normal& n, Attribute& a)
		{
			m_vertices.push_back(v);
			m_normals.push_back(n);
			m_attributes.push_back(a);
		}

		void resizeTriangles(const int size_)
		{
			m_triangles.resize(size_);
		}
		void addTriangle(const int pos, Triangle& t)
		{
			m_triangles[pos] = t;
		}
		void addTriangle(Triangle& t) { m_triangles.push_back(t); }

		void resizeQuadrilaterals(const int size_)
		{
			m_quadrilaterals.resize(size_);
		}

		void addQuadrilateral(const int pos, Quadrilateral& q)
		{
			m_quadrilaterals[pos] = q;
		}

		void addLineVertex(Vertex& v, Normal& n, Attribute& a)
		{
			m_lineVertices.push_back(v);
			m_lineNormals.push_back(n);
			m_lineAttributes.push_back(a);
		}
		void resizeLines(const int size_) {
			m_lines.resize(size_);
		}
		void addLine(Line& l) { m_lines.push_back(l); }
		void addLine(const int pos, const int v0, const int v1) { m_lines[pos] = { v0,v1 }; }

		void resizePoints(const int size_)
		{
			m_points.resize(size_);
			m_pointNormals.resize(size_);
			m_pointAttributes.resize(size_);
		}
		void addPoint(Point& p, Normal& n)
		{
			m_points.push_back(p);
			m_pointNormals.push_back(n);
			m_pointAttributes.push_back(pointColor);
		}
		void addPoint(Point& p, Normal& n, Attribute& a)
		{
			m_points.push_back(p);
			m_pointNormals.push_back(n);
			m_pointAttributes.push_back(a);
		}

		void setBBox(BBox& bbox) {
			m_bbox = BBox{ bbox };
		}


		auto&& getTriangles() { return m_triangles; }
		auto&& getTriangles() const { return m_triangles; }
		auto&& getQuadrilaterals() { return m_quadrilaterals; }
		auto&& getQuadrilaterals() const { return m_quadrilaterals; }
		auto&& getVertices() { return m_vertices; }
		auto&& getVertices() const { return m_vertices; }
		auto&& getNormals() { return m_normals; }
		auto&& getNormals() const { return m_normals; }
		auto&& getAttributes() { return m_attributes; }
		auto&& getAttributes() const { return m_attributes; }

		void setVertices(const std::vector<Vertex>& v) { std::copy(v.begin(), v.end(), std::back_inserter(m_vertices)); }
		void setNormals(const std::vector<Normal>& n) { std::copy(n.begin(), n.end(), std::back_inserter(m_normals)); }
		void setTriangles(const std::vector<Triangle>& t) { std::copy(t.begin(), t.end(), std::back_inserter(m_triangles)); }
		void setAttributes(const std::vector<Attribute>& a) { std::copy(a.begin(), a.end(), std::back_inserter(m_attributes)); }


        auto trianglesData() { return m_triangles.data()->data(); }
        auto quadrilateralsData() { return m_quadrilaterals.data()->data(); }
        auto verticesData() { return m_vertices[0].data(); }
        auto normalsData() { return m_normals[0].data(); }
        auto linesData() { return m_lines.data()->data(); }

		auto&& getLineVertices() { return m_lineVertices; }
		auto&& getLineNormals() { return m_lineVertices; }
		auto&& getLineAttributes() { return m_lineVertices; }
		auto&& getLines() { return m_lines; }

		Vertex vertex(const int i) { return m_vertices[i]; }
		const Vertex vertex(const int i) const { return m_vertices[i]; }
		void vertex(const int i, const Vertex& v) { m_vertices[i] = v; }

		Normal normal(const int i) { return m_normals[i]; }
		const Normal normal(const int i) const { return m_normals[i]; }

		Triangle triangle(const int i) { return m_triangles[i]; }
		const Triangle triangle(const int i) const { return m_triangles[i]; }

		Quadrilateral quadrilateral(const int i) { return m_quadrilaterals[i]; }
		Quadrilateral quadrilateral(const int i) const { return m_quadrilaterals[i]; }

		int quadVertex(const int q, const int i) { return m_quadrilaterals[q][i]; }
		int quadVertex(const int q, const int i) const { return m_quadrilaterals[q][i]; }

		Attribute attribute(const int i) { return m_attributes[i]; }
		const Attribute attribute(const int i) const { return m_attributes[i]; }
		void attribute(const int i, Attribute& a) {
			m_attributes[i] = a;
		}

		Vertex lineVertex(const int i) { return m_lineVertices[i]; }
		Normal lineNormal(const int i) { return m_lineNormals[i]; }
		Attribute lineAttribute(const int i) { return m_lineAttributes[i]; }
		const Vertex lineVertex(const int i) const { return m_lineVertices[i]; }
		const Normal lineNormal(const int i) const { return m_lineNormals[i]; }
		const Attribute lineAttribute(const int i) const { return m_lineAttributes[i]; }
		size_t nrLineVertices() { return m_lineVertices.size(); }
		Line line(const int i) { return m_lines[i]; }
		const Line line(const int i) const { return m_lines[i]; }

		Point point(const int i) { return m_points[i]; }
		const Point point(const int i) const { return m_points[i]; }
		Normal pointNormal(const int i) { return m_pointNormals[i]; }
		const Normal pointNormal(const int i) const { return m_pointNormals[i]; }
		Attribute pointAttribute(const int i) { return m_pointAttributes[i]; }
		const Attribute pointAttribute(const int i) const { return m_pointAttributes[i]; }

		// modify elements
		void attribute(const int index, const Attribute& a) {
			m_attributes[index] = a;
		}

		size_t nrVertices() { return m_vertices.size(); }
		size_t nrVertices() const { return m_vertices.size(); }
		size_t nrTriangles() { return m_triangles.size(); }
		size_t nrTriangles() const { return m_triangles.size(); }
		size_t nrQuadrilaterals() { return m_quadrilaterals.size(); }
		size_t nrQuadrilaterals() const { return m_quadrilaterals.size(); }
		size_t nrLines() { return m_lines.size(); }
		size_t nrLines() const { return m_lines.size(); }
		size_t nrPoints() { return m_points.size(); }
		size_t nrPoints() const { return m_points.size(); }

		void flipNormals()
		{
			for (auto& n : m_normals)
			{
				n.flip();
			}
		}

		double minX() { return m_bbox[0][0]; }
		double minY() { return m_bbox[0][1]; }
		double minZ() { return m_bbox[0][2]; }
		double maxX() { return m_bbox[7][0]; }
		double maxY() { return m_bbox[7][1]; }
		double maxZ() { return m_bbox[7][2]; }
		double minX() const { return m_bbox[0][0]; }
		double minY() const { return m_bbox[0][1]; }
		double minZ() const { return m_bbox[0][2]; }
		double maxX() const { return m_bbox[7][0]; }
		double maxY() const { return m_bbox[7][1]; }
		double maxZ() const { return m_bbox[7][2]; }
		void bbox()
		{
			double xMin = std::numeric_limits<double>::max();
			double xMax = -xMin;
			double yMin = xMin;
			double yMax = -xMin;
			double zMin = xMin;
			double zMax = -xMin;
			for (auto p : m_vertices)
			{
				xMin = (xMin > p[0]) ? p[0] : xMin;
				xMax = (xMax < p[0]) ? p[0] : xMax;
				yMin = (yMin > p[1]) ? p[1] : yMin;
				yMax = (yMax < p[1]) ? p[1] : yMax;
				zMin = (zMin > p[2]) ? p[2] : zMin;
				zMax = (zMax < p[2]) ? p[2] : zMax;
			}
			m_bbox[0] = Vertex{ xMin, yMin, zMin };
			m_bbox[1] = Vertex{ xMax, yMin, zMin };
			m_bbox[2] = Vertex{ xMin, yMax, zMin };
			m_bbox[3] = Vertex{ xMax, yMax, zMin };
			m_bbox[4] = Vertex{ xMin, yMin, zMax };
			m_bbox[5] = Vertex{ xMax, yMin, zMax };
			m_bbox[6] = Vertex{ xMin, yMax, zMax };
			m_bbox[7] = Vertex{ xMax, yMax, zMax };
		}
		void clear()
		{
			m_vertices.clear();
			m_normals.clear();
			m_triangles.clear();
			m_quadrilaterals.clear();
			m_attributes.clear();
			m_lines.clear();
			m_lineVertices.clear();
			m_lineNormals.clear();
			m_lineAttributes.clear();
			m_points.clear();
			m_pointNormals.clear();
			m_pointAttributes.clear();
		}

		void onering()
		{
			if (m_onering.size() != m_vertices.size())
			{
				// the one still have to be computed
				for (auto o : m_onering)
				{
					if (o.empty())
						o.clear();
				}
				m_onering.clear();
				m_onering.resize(m_vertices.size());
				for (auto t : m_triangles)
				{
					m_onering[t[0]].insert(t[1]);
					m_onering[t[0]].insert(t[2]);

					m_onering[t[1]].insert(t[0]);
					m_onering[t[1]].insert(t[2]);

					m_onering[t[2]].insert(t[0]);
					m_onering[t[2]].insert(t[1]);
				}
			}
		}
		OneRing onering(const int v) { return m_onering[v]; }

	private:
		std::vector<Vertex> m_vertices;
		std::vector<Normal> m_normals;
		std::vector<Triangle> m_triangles;
		std::vector<OneRing> m_onering;
		std::vector<Quadrilateral> m_quadrilaterals;
		std::vector<Attribute> m_attributes;

		std::vector<Line> m_lines;
		std::vector<Vertex> m_lineVertices;
		std::vector<Normal> m_lineNormals;
		std::vector<Attribute> m_lineAttributes;

		std::vector<Vertex> m_points;
		std::vector<Normal> m_pointNormals;
		std::vector<Attribute> m_pointAttributes;
		Attribute pointColor{ 72. / 255., 61. / 255., 139. / 255., 1. };

		BBox m_bbox;
	};
}
