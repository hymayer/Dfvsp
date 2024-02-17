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

//�ͷ��ڴ�ռ�
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
            //����һ���µı���
            ArcNode* temp = new ArcNode;
            temp->adjvex = input.adjList[i][j];
            temp->next = NULL;
            //�����ǰ����Ļ�û�б�����ʱ��
            if (this->arc[i].firstarc == NULL) {
                this->arc[i].firstarc = temp;
            }
            else {
                ArcNode* now = this->arc[i].firstarc;
                while (now->next) {
                    now = now->next;
                }//�ҵ�����������һ�����
                now->next = temp;
            }
            ++count;
        }
    }
}

void Graph_DG::print() {
    int count = 0;
    cout << "ͼ���ڽӾ���Ϊ��" << endl;
    //��������������������
    while (count != this->vexnum) {
        //�������Ľ��
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
    //ջs���ڱ���ջΪ�յĶ����±�
    stack<int> s;
    int i;
    ArcNode* temp;
    //����ÿ���������ȣ�������indgree������
    for (i = 0; i != this->vexnum; i++) {
        temp = this->arc[i].firstarc;
        while (temp) {
            ++this->indegree[temp->adjvex];
            temp = temp->next;
        }
    }

    //�����Ϊ0�Ķ�����ջ
    for (i = 0; i != this->vexnum; i++) {
        if (!this->indegree[i]) {
            s.push(i);
        }
    }
    //count���ڼ�������Ķ������
    int count = 0;
    while (!s.empty()) {//���ջΪ�գ������ѭ��
        i = s.top();
        s.pop();//����ջ��Ԫ�أ�����ջ��Ԫ�س�ջ
        cout << this->arc[i].data << " ";//�����������
        temp = this->arc[i].firstarc;
        while (temp) {
            if (!(--this->indegree[temp->adjvex])) {//�����ȼ��ٵ�Ϊ0������ջ
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
    return false;//˵�����ͼ�л�
}
