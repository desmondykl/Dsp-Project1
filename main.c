#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//CHANGE BLOCK SIZE HERE
//#define BLOCK 500
#define BLOCK 100

#define MAXRECORD ((BLOCK/20)-1)
//Order
#define order (int)((((BLOCK+4)/12))-2)

// 1 row/record = 20 byte
typedef struct record {
   char  title[12] ;    // 12 btye
   float rating ;       // 4 btye
   int   vote ;         // 4 btye
} Record;

struct Blocks {
    Record records[MAXRECORD];    // MAXRECORD*20 btye
    struct Blocks * next;         // 8 byte
    int size;                     // 4 btye
};

typedef struct link_list_head{
    struct link_list *head;
    struct link_list *last;
}link_list_head;

typedef struct link_list{
    Record *ptr;
	struct link_list *nextItem;
}link_list;

typedef struct node {
    void ** pointers;               // 8*order btye
	float * keys;                   // 4*(order-1) btye
	struct node * parent;           // 8 btye
	bool is_leaf;                   // 1 btye
	int num_keys;                   // 4 btye
	struct node * next;             // 8 btye
} node;

node * insert_into_parent(node * root, node * left, float key, node * right);
node * insert_into_node_after_splitting(node * root, node * parent,int left_index,float key, node * right);
node * delete_entry(node * root, node * n, float key, void * pointer,int *numOfAccess);
node * deleteIndex(node * root, float key,int *numOfAccess);

node * queue = NULL;
bool verbose_output = false;
node *root = NULL;
int numberR = 0;

int height(node * const root) {
	int h = 0;
	node * c = root;
	while (!c->is_leaf) {
		c = c->pointers[0];
		h++;
	}
	return h;
}

// Create a new node
node * make_node(void) {
	node * newNode;
	newNode = malloc(sizeof(node));
	if (newNode == NULL) {
		printf("Error in creating new node");
		exit(EXIT_FAILURE);
	}
	newNode->keys = malloc((order - 1) * sizeof(float));
	if (newNode->keys == NULL) {
		exit(EXIT_FAILURE);
	}
	newNode->pointers = malloc(order * sizeof(void *));
	if (newNode->pointers == NULL) {
		exit(EXIT_FAILURE);
	}
	newNode->is_leaf = false;
	newNode->num_keys = 0;
	newNode->parent = NULL;
	newNode->next = NULL;
	return newNode;
}

node * make_leaf(void) {
	node * leaf = make_node();
	leaf->is_leaf = true;
	return leaf;
}

// initialize a new B+ tree node
node *start_new_tree(float key, link_list *pointer) {
	node * root = make_leaf();
	root->keys[0] = key;
	root->pointers[0] = pointer;
	root->pointers[order - 1] = NULL;
	root->parent = NULL;
	root->num_keys++;
	return root;
}

/* Traces the path from the root to a leaf, searching
 * by key.  Displays information about the path
 * if the verbose flag is set.
 * Returns the leaf containing the given key.
 */
node * find_leaf(node * const root, float key, bool verbose,int *numOfAccess) {
	if (root == NULL) {
		printf("Empty tree.\n");
		return root;
	}
	int i = 0;
	int Counter = 0;
	node * temp = root;
	int check= 0;
	if(numOfAccess != NULL){
         check = (*numOfAccess);
	}
	while (!temp->is_leaf) {
		if (verbose) {
            if(Counter == 0)
                printf("| Root Node     : [");
			else
                printf("| Internal Node : [");
			for (i = 0; i < temp->num_keys - 1; i++)
				printf("%0.2f ", temp->keys[i]);
			printf("%0.2f]\n", temp->keys[i]);
		}
		i = 0;
		while (i < temp->num_keys) {
			if (key >= temp->keys[i]) i++;
			else break;
		}
		//Key position commented first not sure need to print anot
		//if (verbose)
			//printf("%d \n", i);
        //go down the tree until it reach leaf then break;
        if(numOfAccess != NULL)
            (*numOfAccess)++;
		temp = (node *)temp->pointers[i];
		Counter = Counter + 1;
	}
	//Leaf Node
	if (verbose) {
		printf("| Leaf Node     : [");
		for (i = 0; i < temp->num_keys - 1; i++)
			printf("%0.2f ", temp->keys[i]);
		printf("%0.2f] \n", temp->keys[i]);
	}
	return temp;
}


// Insert record pointer into leaf node
node *insert_into_leaf(node * leaf, float key, link_list * pointer) {

	int i, insertion_point;

	insertion_point = 0;
	while (insertion_point < leaf->num_keys && leaf->keys[insertion_point] < key)
		insertion_point++;

	for (i = leaf->num_keys; i > insertion_point; i--) {
		leaf->keys[i] = leaf->keys[i - 1];
		leaf->pointers[i] = leaf->pointers[i - 1];
	}
	leaf->keys[insertion_point] = key;
	leaf->pointers[insertion_point] = pointer;
	leaf->num_keys++;
	return leaf;
}


node * insert_into_node(node * root, node * n, int index, float key, node * right) {
	int i;

	for (i = n->num_keys; i > index; i--) {
		n->pointers[i + 1] = n->pointers[i];
		n->keys[i] = n->keys[i - 1];
	}
	n->pointers[index + 1] = right;
	n->keys[index] = key;
	n->num_keys++;
	return root;

}

//Utility function used to divide the node if it is full
int cut(int length) {
	if (length % 2 == 0)
		return length/2;
	else
		return length/2 + 1;
}

// When node is full, create a new node and split the key
node * insert_into_node_after_splitting(node * root, node * old_node, int index,
		float key, node * right) {

	int i, j, split;
	node * new_node, * child;
	float * temp_keys;
	float k_prime;
	node ** temp_pointers;

	temp_pointers = malloc((order + 1) * sizeof(node *));
	if (temp_pointers == NULL) {
		exit(EXIT_FAILURE);
	}
	temp_keys = malloc(order * sizeof(float));
	if (temp_keys == NULL) {
		exit(EXIT_FAILURE);
	}

	for (i = 0, j = 0; i < old_node->num_keys + 1; i++, j++) {
		if (j == index + 1) j++;
		temp_pointers[j] = old_node->pointers[i];
	}

	for (i = 0, j = 0; i < old_node->num_keys; i++, j++) {
		if (j == index) j++;
		temp_keys[j] = old_node->keys[i];
	}

	temp_pointers[index + 1] = right;
	temp_keys[index] = key;

	//Create a new node
	split = cut(order);
	new_node = make_node();
	old_node->num_keys = 0;
	for (i = 0; i < split - 1; i++) {
		old_node->pointers[i] = temp_pointers[i];
		old_node->keys[i] = temp_keys[i];
		old_node->num_keys++;
	}
	old_node->pointers[i] = temp_pointers[i];
	k_prime = temp_keys[split - 1];
	for (++i, j = 0; i < order; i++, j++) {
		new_node->pointers[j] = temp_pointers[i];
		new_node->keys[j] = temp_keys[i];
		new_node->num_keys++;
	}
	new_node->pointers[j] = temp_pointers[i];
	free(temp_pointers);
	free(temp_keys);
	new_node->parent = old_node->parent;
	for (i = 0; i <= new_node->num_keys; i++) {
		child = new_node->pointers[i];
		child->parent = new_node;
	}
    //Update the key in the parent node
	return insert_into_parent(root, old_node, k_prime, new_node);
}

node * insert_into_new_root(node * left, float key, node * right) {
	node * root = make_node();
	root->keys[0] = key;
	root->pointers[0] = left;
	root->pointers[1] = right;
	root->num_keys++;
	root->parent = NULL;
	left->parent = root;
	right->parent = root;
	return root;
}

int get_left_index(node * parent, node * left) {

	int left_index = 0;
	while (left_index <= parent->num_keys &&
			parent->pointers[left_index] != left)
		left_index++;
	return left_index;
}

node * insert_into_parent(node * root, node * left, float key, node * right) {
	int index;
	node * parent;
	parent = left->parent;
	if (parent == NULL)
		return insert_into_new_root(left, key, right);
    index = get_left_index(parent, left);
    if (parent->num_keys < order - 1)
		return insert_into_node(root, parent, index, key, right);
	return insert_into_node_after_splitting(root, parent, index, key, right);
}

node * insert_into_leaf_after_splitting(node * root, node * leaf, float key, link_list * pointer) {

	node * new_leaf;
	float * temp_keys;
	void ** temp_pointers;
	int insertion_index, split, i, j;
	float new_key;

	new_leaf = make_leaf();

	temp_keys = malloc(order * sizeof(float));
	if (temp_keys == NULL) {
		exit(EXIT_FAILURE);
	}

	temp_pointers = malloc(order * sizeof(link_list));
	if (temp_pointers == NULL) {
		exit(EXIT_FAILURE);
	}

	insertion_index = 0;
	while (insertion_index < order - 1 && leaf->keys[insertion_index] < key)
		insertion_index++;
	for (i = 0, j = 0; i < leaf->num_keys; i++, j++) {
		if (j == insertion_index) j++;
		temp_keys[j] = leaf->keys[i];
		temp_pointers[j] = leaf->pointers[i];
	}

	temp_keys[insertion_index] = key;
	temp_pointers[insertion_index] = pointer;

	leaf->num_keys = 0;

	split = cut(order - 1);

	for (i = 0; i < split; i++) {
		leaf->pointers[i] = temp_pointers[i];
		leaf->keys[i] = temp_keys[i];
		leaf->num_keys++;
	}

	for (i = split, j = 0; i < order; i++, j++) {
		new_leaf->pointers[j] = temp_pointers[i];
		new_leaf->keys[j] = temp_keys[i];
		new_leaf->num_keys++;
	}

	free(temp_pointers);
	free(temp_keys);

	new_leaf->pointers[order - 1] = leaf->pointers[order - 1];
	leaf->pointers[order - 1] = new_leaf;
	for (i = leaf->num_keys; i < order - 1; i++){
		leaf->pointers[i] = NULL;
    }
	for (i = new_leaf->num_keys; i < order - 1; i++)
		new_leaf->pointers[i] = NULL;

	new_leaf->parent = leaf->parent;
	new_key = new_leaf->keys[0];

	return insert_into_parent(root, leaf, new_key, new_leaf);
}

//check for node height
int path_to_root(node * const root, node * child) {
	int length = 0;
	node * c = child;
	while (c != root) {
		c = c->parent;
		length++;
	}
	return length;
}
//Used for printing B+ tree
void enqueue(node * new_node) {
	node * c;
	if (queue == NULL) {
		queue = new_node;
		queue->next = NULL;
	}
	else {
		c = queue;
		while(c->next != NULL) {
			c = c->next;
		}
		c->next = new_node;
		new_node->next = NULL;
	}
}
//Used for printing B+ tree
node * dequeue(void) {
	node * n = queue;
	queue = queue->next;
	n->next = NULL;
	return n;
}

int CalculateLink(node * const root,int count) {
	if (root == NULL) {
		printf("Empty tree.\n");
		return 0;
	}
	int i;
	node * temp = root;
	while (!temp->is_leaf)
		temp = temp->pointers[0];
	while (true) {
		for (i = 0; i < temp->num_keys; i++) {
            link_list_head *lh = temp->pointers[i];
            link_list *ll = lh->head;
            while(ll!=NULL){
                count = count + 1;
                ll = ll->nextItem;
            }
		}
		if (temp->pointers[order - 1] != NULL) {
			temp = temp->pointers[order - 1];
		}
		else
			break;
	}
	return count;
}


void print_tree(node * const root,bool printOutput, int *numOfNode) {

	node * n = NULL;
	int i = 0;
	int rank = 0;
	int new_rank = 0;

	if (root == NULL) {
		printf("Empty tree.\n");
		return;
	}
	queue = NULL;
	enqueue(root);
    if(printOutput)
        printf("Root Node     : ");
	while(queue != NULL) {
		n = dequeue();
		(*numOfNode)++;
		if (n->parent != NULL && n == n->parent->pointers[0]) {
			new_rank = path_to_root(root, n);
			if (new_rank != rank) {
				rank = new_rank;
				if(printOutput){
                    printf("\n");
                    if(!n->is_leaf)
                        printf("Internal Node : ");
                    else
                        printf("Leaf Node     : ");
				}
			}
		}
		for (i = 0; i < n->num_keys; i++) {
            if(i==0)
                if(printOutput)
                    printf("[");
			if(printOutput)
                printf("%0.2f ", n->keys[i]);
		}
		if (!n->is_leaf)
			for (i = 0; i <= n->num_keys; i++)
				enqueue(n->pointers[i]);
		if(printOutput)
            printf("] ");
	}
}

// Look for record ptr based on the given key
link_list * find(node * root, float key, bool verbose, node ** leaf_out, int *numOfAccess) {
    if (root == NULL) {
        if (leaf_out != NULL) {
            *leaf_out = NULL;
        }
        return NULL;
    }

	int i = 0;
    node * leaf = NULL;
	leaf = find_leaf(root, key, verbose,numOfAccess);

	for (i = 0; i < leaf->num_keys; i++)
		if (leaf->keys[i] == key) break;
    if (leaf_out != NULL) {
        *leaf_out = leaf;
    }
	if (i == leaf->num_keys)
		return NULL;
	else
		return (link_list*)leaf->pointers[i];
}

int find_range(node * const root, float key_start, float key_end, bool verbose, int *numOfAccess,int *numOfHead) {
	int i, num_found,k;
	num_found = 0;
	k=0;
	node * n =NULL;
	if(verbose)
        n = find_leaf(root, key_start, false,numOfAccess);
    else
        n = find_leaf(root, key_start, true,numOfAccess);

	if (n == NULL) return 0;
	//get position of pointer in the leaf node
	for (i = 0; i < n->num_keys && n->keys[i] < key_start; i++);
	if (i == n->num_keys) return 0;

	int g=i;
	node *print = n->pointers[order - 1];
	if(!verbose){
        while (print != NULL) {
            bool toPrint = false;
            for (int g = 0; g < print->num_keys && print->keys[g] <= key_end; g++) {
                toPrint=true;
                break;
            }
            if(toPrint){
                int a= 0;
                printf("| Leaf Node     : [");
                for (a = 0; a < print->num_keys - 1; a++)
                    printf("%0.2f ", print->keys[a]);
                printf("%0.2f]\n", print->keys[a]);
                (*numOfAccess)++;
            }
            print = print->pointers[order - 1];
        }
	}


    if(verbose){
        printf("+-----------------------------------------------------------+\n");
        printf("|        tconst of Record Found (First 10 Records)          |\n");
        printf("+-----------+-----------+-----------+-----------+-----------+\n");
    }
	while (n != NULL) {
		for (; i < n->num_keys && n->keys[i] <= key_end; i++) {
            link_list_head *llhead = n->pointers[i];
            link_list *ll = llhead->head;
            k++;
            while(ll != NULL){
                if(verbose && num_found <10 ){
                    printf("| %s " , ll->ptr->title);
                }
                if(verbose && num_found==10) return;
                ll = ll->nextItem;
                num_found++;
                if(verbose && num_found%5==0)
                    printf("|\n");
            }
		}
		n = n->pointers[order - 1];
		i = 0;
	}
	return num_found;
}

// Insert as a linklist if the key is duplicate
insertAsLink(link_list_head *duplicateKey,Record *ptr){
    link_list_head *ll = duplicateKey;

    link_list *ll_new = (link_list *)malloc(sizeof(link_list));
    ll_new->ptr = ptr;
    ll_new->nextItem=NULL;
    ll->last->nextItem=ll_new;
    ll->last=ll_new;

}


node *insert(node * root,float key,Record *ptr) {
    link_list *record_pointer = NULL;
    node *leaf = NULL;
    link_list *duplicateKey = NULL;

    if (root == NULL){
        link_list *ll = (link_list *)malloc(sizeof(link_list));
        ll->ptr = ptr;
        ll->nextItem=NULL;
        record_pointer=ll;

        link_list_head *llHead = (link_list_head *)malloc(sizeof(link_list_head));
        llHead->head=ll;
        llHead->last=ll;

		return start_new_tree(key, llHead);
    }

    duplicateKey = find(root, key, false, NULL, NULL);

    if(duplicateKey != NULL){
        insertAsLink(duplicateKey, ptr);
        return root;
    }

    leaf = find_leaf(root, key, false,NULL);
    link_list *ll = (link_list *)malloc(sizeof(link_list));
	ll->ptr = ptr;
	ll->nextItem=NULL;
    record_pointer=ll;

    link_list_head *llHead = (link_list_head *)malloc(sizeof(link_list_head));
    llHead->head=ll;
    llHead->last=ll;
    record_pointer = llHead;
    if (leaf->num_keys < order - 1) {
		leaf = insert_into_leaf(leaf, key, record_pointer);
		return root;
	}
	/* Case:  leaf must be split.
	 */
	return insert_into_leaf_after_splitting(root, leaf, key, record_pointer);

}

struct Blocks *insertData(struct Blocks *head, char *line){
    struct Blocks *current = head;
    struct Blocks *prev = NULL;
    while (current != NULL) {
        if(current->size < MAXRECORD)
            break;
        prev = current;
        current = current->next;
    }
    if(current==NULL){
        struct Blocks *newBlock = (struct Blocks *)malloc(sizeof(struct Blocks));
        current=newBlock;
        current->size=0;
        current->next=NULL;
        prev->next=newBlock;
    }
    char *token;
    token = strtok(line, "\t");
    Record *record = malloc(sizeof(Record));

    int type =0;
    float rating_key = 0;
    while( token != NULL ) {
        // each token is a value for the line
        if(type==0)
            strcpy(record->title, token);
        else if(type==1){
            record->rating = atof(token);
            rating_key = atof(token);
        }
        else{
            record->vote = atol(token);
        }
        type++;
        token = strtok(NULL, "\t");
    }
    current->records[current->size] = *record;
    root = insert(root, rating_key,&current->records[current->size]);

    current->size++;
    return current;
}

//pretty print data block
void print2ndBlockRecord (Record record ){
    if(record.vote<=9)
        printf(" %s  | %0.2f   | %d     |\n", record.title,record.rating,record.vote);
    else if(record.vote<=99)
        printf(" %s  | %0.2f   | %d    |\n", record.title,record.rating,record.vote);
    else if(record.vote<=999)
        printf(" %s  | %0.2f   | %d   |\n", record.title,record.rating,record.vote);
    else if(record.vote<=9999)
        printf(" %s  | %0.2f   | %d  |\n", record.title,record.rating,record.vote);
    else if(record.vote<=99999)
        printf(" %s  | %0.2f   | %d |\n", record.title,record.rating,record.vote);
    else if(record.vote<=999999)
        printf(" %s  | %0.2f   | %d|\n", record.title,record.rating,record.vote);
    else
        printf(" %s  | %0.2f   | %d   |\n", record.title,record.rating,record.vote);
}
//pretty print data block
void printDataBlock(struct Blocks *blocks[]){
    struct Blocks *Blocks = blocks[0];
    struct Blocks *Blocks2 = blocks[1];
    for(int i = 0; i<MAXRECORD; i++){

        if(Blocks->records[i].vote<=9){
            printf("| %s  | %0.2f   | %d     |", Blocks->records[i].title,Blocks->records[i].rating,Blocks->records[i].vote);
            print2ndBlockRecord(Blocks2->records[i]);
        }
        else if(Blocks->records[i].vote<=99){
            printf("| %s  | %0.2f   | %d    |", Blocks->records[i].title,Blocks->records[i].rating,Blocks->records[i].vote);
            print2ndBlockRecord(Blocks2->records[i]);
        }
        else if(Blocks->records[i].vote<=999){
            printf("| %s  | %0.2f   | %d   |", Blocks->records[i].title,Blocks->records[i].rating,Blocks->records[i].vote);
            print2ndBlockRecord(Blocks2->records[i]);
        }
        else if(Blocks->records[i].vote<=9999){
            printf("| %s  | %0.2f   | %d  |", Blocks->records[i].title,Blocks->records[i].rating,Blocks->records[i].vote);
            print2ndBlockRecord(Blocks2->records[i]);
        }
        else if(Blocks->records[i].vote<=99999){
            printf("| %s  | %0.2f   | %d |", Blocks->records[i].title,Blocks->records[i].rating,Blocks->records[i].vote);
            print2ndBlockRecord(Blocks2->records[i]);
        }
        else if(Blocks->records[i].vote<=999999){
            printf("| %s  | %0.2f   | %d|", Blocks->records[i].title,Blocks->records[i].rating,Blocks->records[i].vote);
            print2ndBlockRecord(Blocks2->records[i]);
        }
        else{
            printf("| %s  | %0.2f   | %d   |", Blocks->records[i].title,Blocks->records[i].rating,Blocks->records[i].vote);
            print2ndBlockRecord(Blocks2->records[i]);
        }

        if(i==MAXRECORD-1)
            printf("+------------+--------+-------+------------+--------+-------+\n");
    }
}
int calculateblockAccessRange (float sKey,float eKey,struct Blocks *head, bool verbose,int numOfHead){

    int numofaccess=0;
    struct Blocks *currentBlock = head;
    struct Blocks *oldBlock = NULL;
    struct Blocks *printBlock[2] ;
    int p = 0;
    while(currentBlock!=NULL){
        bool sameKey = false;
        float arraykey[MAXRECORD];
        int c = 0;
        bool toggle = true;

        for(int i = 0; i<MAXRECORD; i++){
            Record *currentBlockRecord = &currentBlock->records[i];
            if(currentBlockRecord->rating>= sKey && currentBlockRecord->rating<= eKey){
                for(int k = 0; k<MAXRECORD; k++){
                    if(arraykey[k] == currentBlockRecord->rating){
                        sameKey=true;
                        break;
                    }
                    else{
                        arraykey[c] = currentBlockRecord->rating;
                        c = c+1;
                        numofaccess++;
                        break;
                    }
                }

                if(toggle){
                    printBlock[p] = currentBlock;
                    p++;
                    toggle=false;
                }
            }
        }
        if(p==2){
            if(verbose && numofaccess<11)
                printDataBlock(printBlock);
            p=0;
            printBlock[0] = NULL;
            printBlock[1] = NULL;
        }
        currentBlock = currentBlock->next;
    }
    return numofaccess;

}

int calculateblockAccess (link_list_head *ll,struct Blocks *head,bool verbose,float Key){

    if(ll==NULL) return;

    int numofaccess=0;
    struct Blocks *currentBlock = head;
    int p = 0;
    struct Blocks *printBlock[2] ;
    while(currentBlock!=NULL){
        bool sameKey = false;
        bool toggle = true;
        for(int i = 0; i<MAXRECORD; i++){
            Record *currentBlockRecord = &currentBlock->records[i];
            if(currentBlockRecord->rating == Key){
                if(!sameKey)
                    numofaccess++;
                sameKey = true;
                if(toggle){
                    printBlock[p] = currentBlock;
                    p++;
                    toggle=false;
                }
                //break;
            }
        }
        if(p==2){
            if(verbose && numofaccess<11)
                printDataBlock(printBlock);
            p=0;
            printBlock[0] = NULL;
            printBlock[1] = NULL;
        }
        currentBlock = currentBlock->next;
    }
    return numofaccess;

}

node * remove_entry_from_node(node * n, float key, node * pointer) {

	int i, num_pointers;

	i = 0;
	while (n->keys[i] != key)
		i++;
	for (++i; i < n->num_keys; i++)
		n->keys[i - 1] = n->keys[i];

	num_pointers = n->is_leaf ? n->num_keys : n->num_keys + 1;
	i = 0;
	while (n->pointers[i] != pointer)
		i++;
	for (++i; i < num_pointers; i++)
		n->pointers[i - 1] = n->pointers[i];

	n->num_keys--;
	if (n->is_leaf)
		for (i = n->num_keys; i < order - 1; i++)
			n->pointers[i] = NULL;
	else
		for (i = n->num_keys + 1; i < order; i++)
			n->pointers[i] = NULL;
	return n;
}

node * adjust_root(node * root) {

	node * new_root;
	if (root->num_keys > 0)
		return root;
	if (!root->is_leaf) {
		new_root = root->pointers[0];
		new_root->parent = NULL;
	}
	else
		new_root = NULL;

	free(root->keys);
	free(root->pointers);
	free(root);
	return new_root;
}

int get_neighbor_index(node * n) {
	int i;
	for (i = 0; i <= n->parent->num_keys; i++)
		if (n->parent->pointers[i] == n)
			return i - 1;
	exit(EXIT_FAILURE);
}

node * coalesce_nodes(node * root, node * n, node * neighbor, int neighbor_index, float k_prime,int *numOfaccess) {

	int i, j, neighbor_insertion_index, n_end;
	node * tmp;
	if (neighbor_index == -1) {
		tmp = n;
		n = neighbor;
		neighbor = tmp;
	}
	neighbor_insertion_index = neighbor->num_keys;
	if (!n->is_leaf) {
        neighbor->keys[neighbor_insertion_index] = k_prime;
		neighbor->num_keys++;
        n_end = n->num_keys;
        for (i = neighbor_insertion_index + 1, j = 0; j < n_end; i++, j++) {
			neighbor->keys[i] = n->keys[j];
			neighbor->pointers[i] = n->pointers[j];
			neighbor->num_keys++;
			n->num_keys--;
		}
        neighbor->pointers[i] = n->pointers[j];
        for (i = 0; i < neighbor->num_keys + 1; i++) {
			tmp = (node *)neighbor->pointers[i];
			tmp->parent = neighbor;
		}
	}
	else {
		for (i = neighbor_insertion_index, j = 0; j < n->num_keys; i++, j++) {
			neighbor->keys[i] = n->keys[j];
			neighbor->pointers[i] = n->pointers[j];
			neighbor->num_keys++;
		}
		neighbor->pointers[order - 1] = n->pointers[order - 1];
	}
	root = delete_entry(root, n->parent, k_prime, n,numOfaccess);
	free(n->keys);
	free(n->pointers);
	free(n);
	(*numOfaccess)++;
	return root;
}

node * redistribute_nodes(node * root, node * n, node * neighbor, int neighbor_index,
		int k_prime_index, float k_prime) {

	int i;
	node * tmp;
	if (neighbor_index != -1) {
		if (!n->is_leaf)
			n->pointers[n->num_keys + 1] = n->pointers[n->num_keys];
		for (i = n->num_keys; i > 0; i--) {
			n->keys[i] = n->keys[i - 1];
			n->pointers[i] = n->pointers[i - 1];
		}
		if (!n->is_leaf) {
			n->pointers[0] = neighbor->pointers[neighbor->num_keys];
			tmp = (node *)n->pointers[0];
			tmp->parent = n;
			neighbor->pointers[neighbor->num_keys] = NULL;
			n->keys[0] = k_prime;
			n->parent->keys[k_prime_index] = neighbor->keys[neighbor->num_keys - 1];
		}
		else {
			n->pointers[0] = neighbor->pointers[neighbor->num_keys - 1];
			neighbor->pointers[neighbor->num_keys - 1] = NULL;
			n->keys[0] = neighbor->keys[neighbor->num_keys - 1];
			n->parent->keys[k_prime_index] = n->keys[0];
		}
	}
	else {
		if (n->is_leaf) {
			n->keys[n->num_keys] = neighbor->keys[0];
			n->pointers[n->num_keys] = neighbor->pointers[0];
			n->parent->keys[k_prime_index] = neighbor->keys[1];
		}
		else {
			n->keys[n->num_keys] = k_prime;
			n->pointers[n->num_keys + 1] = neighbor->pointers[0];
			tmp = (node *)n->pointers[n->num_keys + 1];
			tmp->parent = n;
			n->parent->keys[k_prime_index] = neighbor->keys[0];
		}
		for (i = 0; i < neighbor->num_keys - 1; i++) {
			neighbor->keys[i] = neighbor->keys[i + 1];
			neighbor->pointers[i] = neighbor->pointers[i + 1];
		}
		if (!n->is_leaf)
			neighbor->pointers[i] = neighbor->pointers[i + 1];
	}

	n->num_keys++;
	neighbor->num_keys--;

	return root;
}

node * delete_entry(node * root, node * n, float key, void * pointer, int *numOfAccess) {

	float min_keys;
	node * neighbor;
	int neighbor_index;
	int k_prime_index;
	float k_prime;
	int capacity;
	n = remove_entry_from_node(n, key, pointer);

	if (n == root)
		return adjust_root(root);

	min_keys = n->is_leaf ? cut(order - 1) : cut(order) - 1;

	if (n->num_keys >= min_keys)
		return root;

	neighbor_index = get_neighbor_index(n);
	k_prime_index = neighbor_index == -1 ? 0 : neighbor_index;
	k_prime = n->parent->keys[k_prime_index];
	neighbor = neighbor_index == -1 ? n->parent->pointers[1] :
		n->parent->pointers[neighbor_index];

	capacity = n->is_leaf ? order : order - 1;
	if (neighbor->num_keys + n->num_keys < capacity){
		return coalesce_nodes(root, n, neighbor, neighbor_index, k_prime,numOfAccess);
    }
	else
		return redistribute_nodes(root, n, neighbor, neighbor_index, k_prime_index, k_prime);
}


node * deleteIndex(node * root, float key,int *numOfAccess) {
	node * key_leaf = NULL;
	link_list_head * key_record = NULL;
	key_record = find(root, key, false, &key_leaf,NULL);
	if (key_record != NULL && key_leaf != NULL) {
		root = delete_entry(root, key_leaf, key, key_record,numOfAccess);
	}
	return root;
}


int main() {

    char line[256];
    FILE * fp;
    const char s[2] = "\t";
    //LOAD THE FILE
    fp = fopen("\data.tsv", "r");
    if (fp == NULL)
        exit(EXIT_FAILURE);

    int numberOfRecordInBlock = 0;
    struct Blocks *head = (struct Blocks *)malloc(sizeof(struct Blocks));
    head->size=0;
    head->next=NULL;
    struct Blocks *last = head;
    while (fgets(line, sizeof(line), fp)) {
        /* note that fgets don't strip the terminating \n, checking its
           presence would allow to handle lines longer that sizeof(line) */
        if(numberR>0){
            last = insertData(last,line);
            //if(numberR==6)break;
        }
        numberR++;
    }

    //Calculate the DataBlock Size
    struct Blocks *current = head;
    int sizeOfDb = 0;
    int NumberOfBlock = 0;
    while (current != NULL) {
        sizeOfDb = sizeof(*current) + sizeOfDb;
        NumberOfBlock++;
        current = current->next;
    }
    fclose(fp);

    //Calculate Index size
    int numberOfNode = 0;
    print_tree(root,false,&numberOfNode);
	int sizeOfNode= (4*(order-1)) + 8+8+4+1+ (8*order);
	int sizeOfIndex = sizeOfNode*numberOfNode;

    int numberOfLink;
    numberOfLink = CalculateLink(root,numberOfLink);
    int sizeOfLinklist = numberOfLink*16;

    printf("\n          DATABASE STATISTICS     \n");
    printf("+------------------------------------+\n");
    printf("| Size of a Block      : %d Btye    |\n",BLOCK);
    printf("| The Size of Database : %d Mb      |\n",  (sizeOfDb+sizeOfIndex+sizeOfLinklist)/100000 );
    printf("| Number of Blocks     : %d\t     |\n",NumberOfBlock);
    printf("| Number of Records    : %d     |\n",numberR);
    printf("+------------------------------------+\n\n");

    //false = just print stats - numOfnodes
    //true = print whole tree;
    int h = height(root);
    numberOfNode = 0;
    printf("\n       B+ TREE STATISTICS\n");
    printf("+-------------------------------+\n");
    printf("| Parameter n         : %d\t|\n", order-1);
    print_tree(root,false,&numberOfNode);
    printf("| Number Of Nodes     : %d\t|\n", numberOfNode);
    printf("| Height of B+ Tree   : %d\t|\n", h);
    printf("+-------------------------------+\n\n");
    printf("\n                  B+ TREE CONTENT\n");
    printf("+-----------------------------------------------------------+\n");
    print_tree(root,true,&numberOfNode);
    printf("\n+-----------------------------------------------------------+\n");

    //Retrieve for a single key
    //Hard Code searching for key value 8 -> can change to let user type the value
    printf("\n\n               SINGLE KEY SEARCH STATISTICS \n");
    printf("+----------------------------------------------------------+\n");
    int numOfAccess = 0 ;
    int numOfRecordFound = 0 ;
    float Key = 8;
    link_list_head *llhead = find(root, Key, true, NULL,&numOfAccess);
    printf("+-----------------------------------------------------------+\n");
    printf("| Number of index block access : %d                          |\n",numOfAccess+1);
    if(llhead != NULL){
        link_list *ll = llhead->head;
        while(ll!=NULL){
            numOfRecordFound++;
            ll=ll->nextItem;
        }
    }
    else{
        printf("No Record Found ! \n");
    }
    printf("| Number of Record Found       : %d                      |\n",numOfRecordFound);
    int blockaccess = calculateblockAccess(llhead,head,false,Key);
    printf("| Number of data block access  : %d                      |\n",blockaccess);
    printf("+-----------------------------------------------------------+\n");
    printf("|        tconst of Record Found (First 10 Records)          |\n");
    printf("+-----------+-----------+-----------+-----------+-----------+\n");
    link_list *ll = llhead->head;
    numOfRecordFound = 0;
    while(ll!=NULL){
        if(numOfRecordFound<10){
            if(numOfRecordFound<9)
                printf("| %s " , ll->ptr->title);
            else{
                printf("| %s |" , ll->ptr->title);
                break;
            }
        }
        numOfRecordFound++;
        ll=ll->nextItem;
        if(numOfRecordFound%5==0)
            printf("|\n");
    }
    printf("\n");
    printf("+-----------+-----------+-----------+-----------+-----------+\n");
    printf("|           Data Block Content (First 10 Accessed)         |\n");
    printf("+------------+--------+-------+------------+--------+-------+\n");
    calculateblockAccess(llhead,head,true,Key);


    printf("\n\n             RANGE OF KEY SEARCH STATISTICS  \n");
    printf("+--------------------------------------------------------+\n");
    //Retrieve for range of key
    //Hard Code searching for range value -> can change to let user type the value
    numOfAccess = 0 ;
    blockaccess = 0;
    int numOfHead = 0;
	float sKey = 7;
	float eKey = 9;

	numOfRecordFound = find_range(root, sKey, eKey, false,&numOfAccess,&numOfHead);
    printf("+-----------------------------------------------------------+\n");
    printf("| Number of Index Block Access: %d                           |\n",numOfAccess+1);
    printf("| Number of Record Found      : %d                      |\n",numOfRecordFound);
    blockaccess = calculateblockAccessRange(sKey, eKey,head,false,numOfHead);
    printf("| Number of Data Block Access : %d                      |\n",blockaccess);
    find_range(root, sKey, eKey, true,&numOfAccess,&numOfHead);
    printf("+-----------+-----------+-----------+-----------+-----------+\n");
    printf("|           Data Block Content (First 10 Accessed)          |\n");
    printf("+------------+--------+-------+------------+--------+-------+\n");
    printf("|   tconst   | Rating | Votes |   tconst   | Rating | Votes |\n");
    printf("+------------+--------+-------+------------+--------+-------+\n");
    numOfHead = 0;
    //Print Block
    calculateblockAccessRange(sKey, eKey,head,true,numOfHead);
    numOfAccess = 0 ;
    blockaccess = 0;

    printf("\n\n    DELETE RECORD STATISTICS\n");
    printf("+-----------------------------------+\n");
    numOfAccess = 0 ;
    deleteIndex(root,7,&numOfAccess);
    printf("| Number of Index Node Deleted : %d  |\n",numOfAccess);
    h = height(root);
    printf("| Height of B+ Tree            : %d  |\n", h);
    printf("+-----------------------------------+\n");
    printf("\n\n                 B+ TREE AFTER DELETION\n");
    printf("+-----------------------------------------------------------+\n");
    print_tree(root,true,&numOfAccess);
    printf("\n+-----------------------------------------------------------+\n");

    printf("\n\n      --END--\n");
}


