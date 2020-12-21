#ifndef __LINKED_LIST_H__
#define __LINKED_LIST_H__

template <typename T>
class LinkedListNode {
    private:
    LinkedListNode ** prev = NULL;
    LinkedListNode * nxt = NULL;
    public:
    T payload;
    inline LinkedListNode(T _payload):payload(_payload){};
    inline void remove() {
        if (prev) {
            *prev = nxt;
        }
        if (nxt) {
            nxt->prev = prev;
        }
        nxt = NULL;
        prev = NULL;
    }
    inline void insert(LinkedListNode *& pos) {
        this->nxt = pos;
        pos = this;
        this->prev = &pos;
        if (this->nxt) {
            this->nxt->prev = &(this->nxt);
        }
    }
    inline LinkedListNode* next() const{
        return this->nxt;
    }
};

#endif