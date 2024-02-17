#include"graph.h"

using namespace szx;

Graph_DG::Graph_DG(int vexnum, int edge) {
    this->vexnum = vexnum;
    this->edge = edge;
    this->arc = new Vnode[this->vexnum];
    this->indegree = new int[this->vexnum];
    for (int i = 0; i < this->vexnum; i++) {
        this->indegree[i] = 0;
        this->arc[i].firstarc = NULL;
        this->arc[i].data = "v" + to_string(i);
    }
}

//释放内存空间
Graph_DG::~Graph_DG() {
    ArcNode* p, * q;
    for (int i = 0; i < this->vexnum; i++) {
        if (this->arc[i].firstarc) {
            p = this->arc[i].firstarc;
            while (p) {
                q = p->next;
                delete p;
                p = q;
            }
        }
    }
    delete[] this->arc;
    delete[] this->indegree;
}

void Graph_DG::createGraph(DFeedbackVertexSet& input) {
    int count = 0;
    for (int i = 0; i < input.adjList.size(); i++) {
        for (int j = 0; j < input.adjList[i].size(); j++) {
            //声明一个新的表结点
            ArcNode* temp = new ArcNode;
            temp->adjvex = input.adjList[i][j];
            temp->next = NULL;
            //如果当前顶点的还没有边依附时，
            if (this->arc[i].firstarc == NULL) {
                this->arc[i].firstarc = temp;
            }
            else {
                ArcNode* now = this->arc[i].firstarc;
                while (now->next) {
                    now = now->next;
                }//找到该链表的最后一个结点
                now->next = temp;
            }
            ++count;
        }
    }
}

void Graph_DG::print() {
    int count = 0;
    cout << "图的邻接矩阵为：" << endl;
    //遍历链表，输出链表的内容
    while (count != this->vexnum) {
        //输出链表的结点
        cout << this->arc[count].data << " ";
        ArcNode* temp = this->arc[count].firstarc;
        while (temp) {
            cout << "<" << this->arc[count].data << "," << this->arc[temp->adjvex].data << "> ";
            temp = temp->next;
        }
        cout << "^" << endl;
        ++count;
    }
}

bool Graph_DG::topological_sort() {
    cout << "Topological sort is:" << endl;
    //栈s用于保存栈为空的顶点下标
    stack<int> s;
    int i;
    ArcNode* temp;
    //计算每个顶点的入度，保存在indgree数组中
    for (i = 0; i != this->vexnum; i++) {
        temp = this->arc[i].firstarc;
        while (temp) {
            ++this->indegree[temp->adjvex];
            temp = temp->next;
        }
    }

    //把入度为0的顶点入栈
    for (i = 0; i != this->vexnum; i++) {
        if (!this->indegree[i]) {
            s.push(i);
        }
    }
    //count用于计算输出的顶点个数
    int count = 0;
    while (!s.empty()) {//如果栈为空，则结束循环
        i = s.top();
        s.pop();//保存栈顶元素，并且栈顶元素出栈
        cout << this->arc[i].data << " ";//输出拓扑序列
        temp = this->arc[i].firstarc;
        while (temp) {
            if (!(--this->indegree[temp->adjvex])) {//如果入度减少到为0，则入栈
                s.push(temp->adjvex);
            }
            temp = temp->next;
        }
        ++count;
    }
    if (count == this->vexnum) {
        cout << endl;
        return true;
    }
    cout << "there is a cycle, no topological sort!" << endl;
    return false;//说明这个图有环
}
