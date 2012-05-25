#include <stdio.h>
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>

using namespace std;

typedef struct _point {
	int x, y;
} point;

typedef struct _line {
	point p0, p1;
} line;

vector<vector<point> > g_buildings; // buildings
vector<point> g_points; // all points
vector<int> g_dist; // distance of each point
vector<int> g_prev; // father id of each point
vector<vector<int> > g_dist_between;
vector<int> g_lineIndex;
int g_n_points;
const int infinity = 0x7FFFFFFF;

int distanceUseFloat(const vector<point> &points)
{
	float f =0;
	point p0, p1;
	p0.x = points[0].x; p0.y = points[0].y;
	for (int i=1; i<points.size(); i++) {
		p1.x = points[i].x; p1.y = points[i].y;
		f += sqrt((float)((p0.x-p1.x)*(p0.x-p1.x)+(p0.y-p1.y)*(p0.y-p1.y)));
		p0 = p1;
	}

	return (int)f;
}

int distanceP(point p0, point p1)
{
	return (int)sqrt((float)((p0.x-p1.x)*(p0.x-p1.x)+(p0.y-p1.y)*(p0.y-p1.y)));
}

// (p1-p0) x (p2-p0)
int xMult(point p0, point p1, point p2)
{
	return (p1.x-p0.x)*(p2.y-p0.y)-(p2.x-p0.x)*(p1.y-p0.y);
}

bool cross(line l0, line l1)
{
	int n0 = xMult(l0.p0, l0.p1, l1.p0);
	int n1 = xMult(l0.p0, l0.p1, l1.p1);
	if (n0*n1 > 0) // ==0 : have same point, will treat as cross.
		return false;

	n0 = xMult(l1.p0, l1.p1, l0.p0);
	n1 = xMult(l1.p0, l1.p1, l0.p1);
	if (n0*n1 > 0)
		return false;
	return true;
}

bool crossLines(line l0, vector<line> lines)
{
	if (lines.empty()) return false;

	for (int i=0; i<(int)lines.size(); i++) {
		if (cross(l0, lines[i])) {
			return true;
		}
	}
	return false;
}

void genCloseLines(vector<line> &lines, const vector<point> points)
{
	if (!lines.empty()) lines.clear();
	line l; 	
	l.p0 = points.back();
	for (int i=0; i<points.size(); i++) {
		l.p1 = points[i];
		lines.push_back(l);
		l.p0 = l.p1;
	}
}

void genRandomLines(vector<line> &lines, const vector<point> points0, const vector<point> points1)
{
	if (!lines.empty()) lines.clear();
	line l;
	if (!points0.empty()) {
		for (int i=0; i<points0.size()-1; i++) {
			l.p0 = points0[i];
			l.p1 = points0[i+1];
			lines.push_back(l);
		}
	}

	if (!points1.empty()) {
		for (int i=0; i<(int)points1.size()-1; i++) {
			l.p0 = points1[i];
			l.p1 = points1[i+1];
			lines.push_back(l);
		}
	}
}

void genPoints(vector<point> &gPoints, const int rm_id, const vector<point> points)
{
	if (!gPoints.empty()) gPoints.clear();

	int ind;
	for (int i=0; i<points.size()-1; i++) {
		ind = (rm_id+1+i)%points.size();
		gPoints.push_back(points.at(ind));
	}
}

void processPoint(const int id, const vector<vector<line> > &lineMap)
{
	// process point
	vector<point> points0, points1;
	line tmpL;
	tmpL.p0 = g_points[id]; // end point
	points1.clear();
	vector<point> bi;
	vector<line> lines;
	int ind_i;
	for (int i=0; i<g_buildings.size(); i++)
	{
		bi = g_buildings.at(i);
		// make a line from point to building.
		for (int i0=0; i0<bi.size(); i0++) {
			tmpL.p1 = bi[i0];
			ind_i = g_lineIndex[i]+i0;
			// check if this line cross the building.
			genPoints(points0, i0, bi);
			genRandomLines(lines, points0, points1);
			if (crossLines(tmpL, lines)) {
				continue;
			}
			// check if this line cross other buildings.
			bool bCross = false;
			for (int k = 0; k<g_buildings.size(); k++) {
				if (k==i) continue;
				lines = lineMap[k];
				if (crossLines(tmpL, lines)) {
					bCross = true;
					break;
				}
			}
			if (!bCross) {
				int tmpDist = distanceP(tmpL.p0, tmpL.p1);
				g_dist_between[id][ind_i] = tmpDist;
				g_dist_between[ind_i][id] = tmpDist;
			}
		}
	}
}

// g_points: 0..n-1 -- building; n -- end; n+1 -- start0; n+2 -- start1; 
void generateGraph()
{
	vector<int> v;
	for (int i=0; i<g_n_points; i++) {
		v.push_back(-1);
	}
	for (int j=0; j<g_n_points; j++) {
		g_dist_between.push_back(v); // -1 -- not visted
	}

	for (int j=0; j<g_n_points; j++) { // 0 -- not connected
		g_dist_between[j][j] = 0;
		g_dist_between[j][g_n_points-3] = 0;
		g_dist_between[j][g_n_points-2] = 0;
		g_dist_between[j][g_n_points-1] = 0;
		g_dist_between[g_n_points-3][j] = 0;
		g_dist_between[g_n_points-2][j] = 0;
		g_dist_between[g_n_points-1][j] = 0;
	}

	// generate line map of each building.
	vector<point> build;
	vector<line> lines;
	vector<vector<line> > lineMap;
	for (int i=0; i<(int)g_buildings.size(); i++) {
		build = g_buildings.at(i);
		genCloseLines(lines, build);
		lineMap.push_back(lines);
	}

	// process points between buildings
	line tmpL;
	int s0, s1;
	s0= s1 = 0;
	bool bCross = false;
	int tmpDist = 0;
	vector<point> bi;
	vector<point> bj;
	vector<point> points0, points1;
	int ind_i, ind_j;
	for (int i=0; i<g_buildings.size(); i++) {
		bi = g_buildings.at(i);
		for (int i0=0; i0<bi.size(); i0++) {
			tmpL.p0 = bi[i0];
			for (int j=0; j<g_buildings.size(); j++) {
				bj = g_buildings.at(j);
				// make a line
				for (int j0=0; j0<bj.size(); j0++) {
					tmpL.p1 = bj[j0];
					ind_i = g_lineIndex[i]+i0;
					ind_j = g_lineIndex[j]+j0;
					if (g_dist_between[ind_i][ind_j]>=0) { // this line is processed
						continue;
					}
					g_dist_between[ind_i][ind_j]=0; // init to not-connect
					g_dist_between[ind_j][ind_i]=0; // init to not-connect

					if (i==j) { // same building
						if (abs(i0-j0)>1 && abs(i0-j0)<bi.size()-1)  continue;
					}
					else {
					// check if this line cross the two buildings.
						genPoints(points0, i0, bi);
						genPoints(points1, j0, bj);
						genRandomLines(lines, points0, points1);
						if (crossLines(tmpL, lines)) {
							continue;
						}
					}
					// check if this line cross other buildings.
					bCross = false;
					for (int k = 0; k<g_buildings.size(); k++) {
						if ((k==i)||(k==j)) continue;
						lines = lineMap[k];
						if (crossLines(tmpL, lines)) {
							bCross = true;
							break;
						}
					}
					if (!bCross) {
						tmpDist = distanceP(tmpL.p0, tmpL.p1);
						g_dist_between[ind_i][ind_j] = tmpDist;
						g_dist_between[ind_j][ind_i] = tmpDist;
					}
				}
			}
		}
	}

	// process point: end, start0, start1
	processPoint(g_n_points-3, lineMap);
	processPoint(g_n_points-2, lineMap);
	processPoint(g_n_points-1, lineMap);
}

bool findPathDijkstra(int end, int start0, int start1)
{
	vector<int> visited;
	for (int i=0; i<g_n_points; i++) {
		g_dist.push_back(infinity);
		g_prev.push_back(-1);
		visited.push_back(0);
	}

	g_dist[end] = 0;
	int current = end;
	int i_shortest = 0; // index of shortest node in neighbor
	int l_shortest;
	bool bFound0 = false;
	bool bFound1 = false;
	int nCount = g_n_points;
	while (nCount--) {
		visited[current] = true;
		l_shortest = infinity;
		//find neighbor
		for (int i=0; i<g_n_points; i++) {
			if (!visited[i] && g_dist_between[current][i]) {
				int dist = g_dist[current]+g_dist_between[current][i];
				if (dist < g_dist[i]) {
					g_dist[i] = dist;
					g_prev[i] = current;
				}
			}
			if (!visited[i] && g_dist[i] < l_shortest) {
				l_shortest=g_dist[i];
				i_shortest=i;
			}
		}
		current = i_shortest;
		if (current==start0) bFound0= true;
		if (current==start1) bFound1= true;
		if (bFound0 && bFound1) return true;
	}
	return false;
}

int main(int argc, char** argv) 
{
//// g_points: 0..n-1 -- building; n -- end; n+1 -- start0; n+2 -- start1; 
	g_n_points = 0;
	g_lineIndex.push_back(0);
	point p;
	vector<point> points;
	ifstream fin("Lost in JingAn Temple.in"); // "test.in"
	fin>>noskipws;
	point p_end, p_s0, p_s1;
	bool bEnd, bS0, bS1;
	bEnd=bS0=bS1=false;
	char ch;
	char str[50];
	while(fin>>ch)
	{
		switch(ch) {
			case '(':
				fin>>p.x;
				fin>>ch;
				fin>>ch;
				fin>>p.y;
				if (!bEnd) {
					bEnd=true;
					p_end = p;
				}
				else if (!bS0) {
					bS0 = true;
					p_s0 = p;
				}
				else if (!bS1) {
					bS1 = true;
					p_s1 = p;
				}
				else {
					points.push_back(p);
					g_points.push_back(p);
				}
				g_n_points++;
				break;
			case '\n':
				if (!points.empty()) {
					g_buildings.push_back(points);
					g_lineIndex.push_back(g_n_points-3); // record index of each line; used when caculate index in g_dist_between.
					points.clear();
				}
				break;
			default: break;
		}
	}
	g_buildings.push_back(points); // push the last line of buildings
	g_points.push_back(p_end);
	g_points.push_back(p_s0);
	g_points.push_back(p_s1);
	generateGraph();
	findPathDijkstra(g_n_points-3,g_n_points-2,g_n_points-1);

	///////////////////////////////////////////////////////////////////////////////////
	// g_dist is not accurate, recaculate with float. 
	vector<point> test;
	test.push_back(p_s0);
	int ind = g_n_points-2;
	while (g_prev[ind] != -1) {
		ind = g_prev[ind];
		//fout<<"(";
		//fout<<g_points[ind].x;
		//fout<<","; fout<<g_points[ind].y;fout<<") ";
		test.push_back(g_points[ind]);
	}

	int s0 = distanceUseFloat(test);
	ind = g_n_points-1;
	test.clear();
	test.push_back(p_s1);
	while (g_prev[ind] != -1) {
		ind = g_prev[ind];
		//fout<<"(";
		//fout<<g_points[ind].x;
		//fout<<","; fout<<g_points[ind].y;fout<<") ";
		test.push_back(g_points[ind]);
	}
	int s1 = distanceUseFloat(test);

	ofstream fout("test.out");
	fout<<"Bus Stop - West Beijing Road, "<<"("<<p_s0.x<<", "<<p_s0.y<<")" <<"\n";
	fout << s0 <<"\n";
	fout<<"Metro YuYuan Road, "<<"("<<p_s1.x<<", "<<p_s1.y<<")" <<"\n";
	fout << s1 <<"\n";

	fin.close();
	fout.close();

	return 0;
}