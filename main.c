// Searching on a B+ Tree in C

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Default order
#define ORDER 3

// 1 row/record = 20 byte
typedef struct record {
   float rating ;
   int   vote ;
   char  title[12] ;
} Record;
//Max 100
struct Blocks {
    Record records[4];
    struct Blocks * next;
    int size;
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
    void ** pointers;
	//void ** pointers;
	float * keys;
	struct node * parent;
	bool is_leaf;
	int num_keys;
	struct node * next; // Used for queue.
} node;

node * insert_into_parent(node * root, node * left, float key, node * right);
node * insert_into_node_after_splitting(node * root, node * parent,int left_index,float key, node * right);
int order = 8;
node * queue = NULL;
bool verbose_output = false;
node *root = NULL;
int numberR = 0;

// NOT NEEDED
Record *make_record(char *title,int vote,float rating) {
	Record * new_record = (Record *)malloc(sizeof(Record));
	if (new_record == NULL) {
		perror("Record creation.");
		exit(EXIT_FAILURE);
	}
	else {
		strcpy(new_record->title, title);
        new_record->rating = rating ;
        new_record->vote = vote;
	}
	return new_record;
}

int height(node * const root) {
	int h = 0;
	node * c = root;
	while (!c->is_leaf) {
		c = c->pointers[0];
		h++;
	}
	return h;
}

/* Creates a new general node, which can be adapted
 * to serve as either a leaf or an internal node.
 */
node * make_node(void) {
	node * new_node;
	new_node = malloc(sizeof(node));
	if (new_node == NULL) {
		perror("Node creation.");
		exit(EXIT_FAILURE);
	}
	new_node->keys = malloc((order - 1) * sizeof(float));
	if (new_node->keys == NULL) {
		perror("New node keys array.");
		exit(EXIT_FAILURE);
	}
	new_node->pointers = malloc(order * sizeof(void *));
	if (new_node->pointers == NULL) {
		perror("New node pointers array.");
		exit(EXIT_FAILURE);
	}
	new_node->is_leaf = false;
	new_node->num_keys = 0;
	new_node->parent = NULL;
	new_node->next = NULL;
	return new_node;
}

node * make_leaf(void) {
	node * leaf = make_node();
	leaf->is_leaf = true;
	return leaf;
}

/* First insertion:
 * start a new tree.
 */
node *start_new_tree(float key, link_list *pointer) {
	node * root = make_leaf();
	root->keys[0] = key;
	root->pointers[0] = pointer;
	root->pointers[order - 1] = NULL;
	root->parent = NULL;
	root->num_keys++;
	//root->nextItem=NULL;
	return root;
}

/* Traces the path from the root to a leaf, searching
 * by key.  Displays information about the path
 * if the verbose flag is set.
 * Returns the leaf containing the given key.
 */
node * find_leaf(node * const root, float key, bool verbose,int *numOfAccess) {
	if (root == NULL) {
		if (verbose)
			printf("Empty tree.\n");
		return root;
	}
	int i = 0;
	node * c = root;
	while (!c->is_leaf) {
		if (verbose) {
			printf("[");
			for (i = 0; i < c->num_keys - 1; i++)
				printf("%f ", c->keys[i]);
			printf("%f]\n", c->keys[i]);
		}
		i = 0;
		while (i < c->num_keys) {
			if (key >= c->keys[i]) i++;
			else break;
		}
		//Key position commented first not sure need to print anot
		//if (verbose)
			//printf("%d \n", i);
        //go down the tree until it reach leaf then break;
        if(numOfAccess != NULL)
            (*numOfAccess)++;
		c = (node *)c->pointers[i];
	}
	//Leaf Node
	if (verbose) {
		printf("Leaf [");
		for (i = 0; i < c->num_keys - 1; i++)
			printf("%f ", c->keys[i]);
		printf("%f] \n", c->keys[i]);
	}
	return c;
}

/* Inserts a new pointer to a record and its corresponding
 * key into a leaf.
 * Returns the altered leaf.
 */
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
	//link_list *ll = (link_list *)malloc(sizeof(link_list));
	//ll->ptr = pointer;
	//ll->nextItem=NULL;
	leaf->pointers[insertion_point] = pointer;
	leaf->num_keys++;
	return leaf;
}
/* Inserts a new key and pointer to a node
 * into a node into which these can fit
 * without violating the B+ tree properties.
 */
node * insert_into_node(node * root, node * n,
		int left_index, float key, node * right) {
	int i;

	for (i = n->num_keys; i > left_index; i--) {
		n->pointers[i + 1] = n->pointers[i];
		n->keys[i] = n->keys[i - 1];
	}
	n->pointers[left_index + 1] = right;
	n->keys[left_index] = key;
	n->num_keys++;
	return root;

}
/* Finds the appropriate place to
 * split a node that is too big into two.
 */
int cut(int length) {
	if (length % 2 == 0)
		return length/2;
	else
		return length/2 + 1;
}
/* Inserts a new key and pointer to a node
 * into a node, causing the node's size to exceed
 * the order, and causing the node to split into two.
 */
node * insert_into_node_after_splitting(node * root, node * old_node, int left_index,
		float key, node * right) {

	int i, j, split;
	node * new_node, * child;
	float * temp_keys;
	float k_prime;
	node ** temp_pointers;

	/* First create a temporary set of keys and pointers
	 * to hold everything in order, including
	 * the new key and pointer, inserted in their
	 * correct places.
	 * Then create a new node and copy half of the
	 * keys and pointers to the old node and
	 * the other half to the new.
	 */

	temp_pointers = malloc((order + 1) * sizeof(node *));
	if (temp_pointers == NULL) {
		perror("Temporary pointers array for splitting nodes.");
		exit(EXIT_FAILURE);
	}
	temp_keys = malloc(order * sizeof(float));
	if (temp_keys == NULL) {
		perror("Temporary keys array for splitting nodes.");
		exit(EXIT_FAILURE);
	}

	for (i = 0, j = 0; i < old_node->num_keys + 1; i++, j++) {
		if (j == left_index + 1) j++;
		temp_pointers[j] = old_node->pointers[i];
	}

	for (i = 0, j = 0; i < old_node->num_keys; i++, j++) {
		if (j == left_index) j++;
		temp_keys[j] = old_node->keys[i];
	}

	temp_pointers[left_index + 1] = right;
	temp_keys[left_index] = key;

	/* Create the new node and copy
	 * half the keys and pointers to the
	 * old and half to the new.
	 */
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

	/* Insert a new key into the parent of the two
	 * nodes resulting from the split, with
	 * the old node to the left and the new to the right.
	 */

	return insert_into_parent(root, old_node, k_prime, new_node);
}
/* Creates a new root for two subtrees
 * and inserts the appropriate key into
 * the new root.
 */
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
/* Helper function used in insert_into_parent
 * to find the index of the parent's pointer to
 * the node to the left of the key to be inserted.
 */
int get_left_index(node * parent, node * left) {

	int left_index = 0;
	while (left_index <= parent->num_keys &&
			parent->pointers[left_index] != left)
		left_index++;
	return left_index;
}
/* Inserts a new node (leaf or internal node) into the B+ tree.
 * Returns the root of the tree after insertion.
 */
node * insert_into_parent(node * root, node * left, float key, node * right) {

	int left_index;
	node * parent;
	parent = left->parent;

	/* Case: new root. */

	if (parent == NULL)
		return insert_into_new_root(left, key, right);

	/* Case: leaf or node. (Remainder of
	 * function body.)
	 */

	/* Find the parent's pointer to the left
	 * node.
	 */

	left_index = get_left_index(parent, left);


	/* Simple case: the new key fits into the node.
	 */

	if (parent->num_keys < order - 1)
		return insert_into_node(root, parent, left_index, key, right);

	/* Harder case:  split a node in order
	 * to preserve the B+ tree properties.
	 */

	return insert_into_node_after_splitting(root, parent, left_index, key, right);
}





/* Inserts a new key and pointer
 * to a new record into a leaf so as to exceed
 * the tree's order, causing the leaf to be split
 * in half.
 */
node * insert_into_leaf_after_splitting(node * root, node * leaf, float key, link_list * pointer) {

	node * new_leaf;
	float * temp_keys;
	void ** temp_pointers;
	int insertion_index, split, i, j;
	float new_key;

	new_leaf = make_leaf();

	temp_keys = malloc(order * sizeof(float));
	if (temp_keys == NULL) {
		perror("Temporary keys array.");
		exit(EXIT_FAILURE);
	}

	temp_pointers = malloc(order * sizeof(link_list));
	if (temp_pointers == NULL) {
		perror("Temporary pointers array.");
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


/* Utility function to give the length in edges
 * of the path from any node to the root.
 */
int path_to_root(node * const root, node * child) {
	int length = 0;
	node * c = child;
	while (c != root) {
		c = c->parent;
		length++;
	}
	return length;
}
/* Helper function for printing the
 * tree out.  See print_tree.
 */
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
/* Helper function for printing the
 * tree out.  See print_tree.
 */
node * dequeue(void) {
	node * n = queue;
	queue = queue->next;
	n->next = NULL;
	return n;
}
/* Prints the bottom row of keys
 * of the tree (with their respective
 * pointers, if the verbose_output flag is set.
 */
void print_leaves(node * const root) {
	if (root == NULL) {
		printf("Empty tree.\n");
		return;
	}
	int i;
	node * c = root;
	while (!c->is_leaf)
		c = c->pointers[0];
	while (true) {
		for (i = 0; i < c->num_keys; i++) {
			if (verbose_output)
				printf("%p ", c->pointers[i]);
			printf("%d ", c->keys[i]);
		}
		if (verbose_output)
			printf("%p ", c->pointers[order - 1]);
		if (c->pointers[order - 1] != NULL) {
			printf(" | ");
			c = c->pointers[order - 1];
		}
		else
			break;
	}
	printf("\n");
}

void print_tree(node * const root,bool printOutput) {

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
	int numberOfNode = 0;

	while(queue != NULL) {
		n = dequeue();
		numberOfNode++;
		if (n->parent != NULL && n == n->parent->pointers[0]) {
			new_rank = path_to_root(root, n);
			if (new_rank != rank) {
				rank = new_rank;
				if(printOutput)
                    printf("\n");
			}
		}
		if (verbose_output)
			printf("(%p)", n);
		for (i = 0; i < n->num_keys; i++) {
			if (verbose_output)
				printf("%p ", n->pointers[i]);
			if(printOutput)
                printf("%0.1f ", n->keys[i]);
		}
		if (!n->is_leaf)
			for (i = 0; i <= n->num_keys; i++)
				enqueue(n->pointers[i]);
		if (verbose_output) {
			if (n->is_leaf)
				printf("%p ", n->pointers[order - 1]);
			else
				printf("%p ", n->pointers[n->num_keys]);
		}
		if(printOutput)
            printf("| ");
	}
	printf("\nNumber Of Nodes : %d", numberOfNode);
	printf("\n");

}

/* Finds and returns the record to which
 * a key refers.
 */
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

    /* If root != NULL, leaf must have a value, even
     * if it does not contain the desired key.
     * (The leaf holds the range of keys that would
     * include the desired key.)
     */

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

/* Finds keys and their pointers, if present, in the range specified
 * by key_start and key_end, inclusive.  Places these in the arrays
 * returned_keys and returned_pointers, and returns the number of
 * entries found.
 */
int find_range(node * const root, float key_start, float key_end, bool verbose,
                 float returned_keys[], void * returned_pointers[], int *numOfAccess,int *numOfHead) {
	int i, num_found,k;
	num_found = 0;
	k=0;
	node * n = find_leaf(root, key_start, true,numOfAccess);
	if (n == NULL) return 0;
	//get position of pointer in the leaf node
	for (i = 0; i < n->num_keys && n->keys[i] < key_start; i++);
	if (i == n->num_keys) return 0;

	int g=i;
	node *print = n->pointers[order - 1];
	while (print != NULL) {
	    bool toPrint = false;
        for (int g = 0; g < print->num_keys && print->keys[g] <= key_end; g++) {
            toPrint=true;
            break;
        }
        if(toPrint){
            int a= 0;
            printf("Leaf [");
            for (a = 0; a < print->num_keys - 1; a++)
                printf("%f ", print->keys[a]);
            printf("%f] \n", print->keys[a]);
        }
        print = print->pointers[order - 1];

	}


	printf("tconst of record found : \n");
	while (n != NULL) {
        (*numOfAccess)++;
		for (; i < n->num_keys && n->keys[i] <= key_end; i++) {
			//returned_keys[num_found] = n->keys[i];
            link_list_head *llhead = n->pointers[i];
            link_list *ll = llhead->head;
            returned_pointers[k] = llhead;
            k++;
            while(ll != NULL){
                if(verbose)
                    printf("%s|" , ll->ptr->title);
                ll = ll->nextItem;
                num_found++;
                //returned_keys[num_found] = n->keys[i];
                //returned_pointers[num_found] = ll->ptr;
            }
		}
		n = n->pointers[order - 1];
		i = 0;
	}
	(*numOfHead)=k;

	return num_found;
}

/* Finds and prints the keys, pointers, and values within a range
 * of keys between key_start and key_end, including both bounds.
 */


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
        if(current->size < 4)
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
void printDataBlock(struct Blocks *Blocks){
    for(int i = 0; i<4; i++){
        if(i==0)
            printf("+----------------+--------+--------+\n");
        if(Blocks->records[i].vote<999)
            printf("| %s\t | %0.2f\t  | %d\t   |\n", Blocks->records[i].title,Blocks->records[i].rating,Blocks->records[i].vote);
        else
            printf("| %s\t | %0.2f\t  | %d   |\n", Blocks->records[i].title,Blocks->records[i].rating,Blocks->records[i].vote);

        if(i==3)
            printf("+----------------+--------+--------+\n");
    }
    printf("\n");
}
int calculateblockAccessRange (float sKey,float eKey,struct Blocks *head, bool verbose,int numOfHead){
    if(verbose)
        printf("Data Block Content :\n");

    int numofaccess=0;
    struct Blocks *currentBlock = head;
    while(currentBlock!=NULL){
        for(int i = 0; i<4; i++){
            Record *currentBlockRecord = &currentBlock->records[i];
            if(currentBlockRecord->rating>= sKey && currentBlockRecord->rating<= eKey){
                numofaccess++;
                if(verbose)
                    printDataBlock(currentBlockRecord);
                break;
            }
        }
        currentBlock = currentBlock->next;
    }
    return numofaccess;

}

int calculateblockAccess (link_list_head *ll,struct Blocks *head,bool verbose){

    if(ll==NULL) return;
    if(verbose)
        printf("Data Block Content :\n");
    //printf("%d",sizeof(link_list));

    struct Blocks *current = head;
    Record *RecordCurrent = ll->head->ptr;
    Record *RecordLast = ll->last->ptr;
    link_list *headLL = ll->head;
    bool found = false;
    int i = 0;

    while (current != NULL && !found) {
        for(i = 0; i<4; i++){
            Record *t = &current->records[i];
            if(t==RecordCurrent){
                found = true;
                break;
            }
        }
        if(found)
            break;
        current = current->next;
    }


    int blockaccess = 0;
    found = false;
    while (current != NULL && !found) {
        for(i = 0; i<4; i++){
            Record *t = &current->records[i];
            if(t == ll->last->ptr){
                if(verbose)
                    printDataBlock(current);
                blockaccess++;
                found = true;
                break;
            }
        }
        if(found)
            break;

        bool match = false;
        for(i = 0; i<4; i++){
            Record *t = &current->records[i];
            if(t==RecordCurrent){
                if(!match){
                    if(verbose)
                        printDataBlock(current);
                    blockaccess++;
                    match = true;
                }
                headLL = headLL->nextItem;
                RecordCurrent = headLL->ptr;
            }


        }

        current = current->next;
    }

    return blockaccess;



}
int main() {


    struct Blocks block;
    printf("Size of double data type : %d\n",sizeof(block));


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

    struct Blocks *current = head;
    int sizeOfDb = 0;
    int NumberOfBlock = 0;
    while (current != NULL) {
        sizeOfDb = sizeof(*current) + sizeOfDb;
        NumberOfBlock++;
        current = current->next;
    }
    //float array[1070319];

    //Never include the index size not sure if need
    printf("\nDatabase statistics\n");
    printf("+------------------------------------+\n");
    printf("| The size of database = %d      |\n",sizeOfDb);
    printf("| Number of Blocks     = %d        |\n",NumberOfBlock);
    printf("| Number of records    = %d       |\n",numberR);
    printf("+------------------------------------+\n\n");

    //false = just print stats - numOfnodes
    //true = print whole tree;
    int h = height(root);
    printf("\nB+ Tree Statistics\n");
    printf("---------------------\n");
    printf("Parameter n = %d ", order);
    print_tree(root,false);
    printf("Height of B+ Tree : %d\n", h);
    printf("---------------------\n\n");


    //Retrieve for a single key
    //Hard Code searching for key value 8 -> can change to let user type the value
    printf("\nSingle Key Search Statistics\n");
    printf("-----------------------------------\n");
    int numOfAccess = 0 ;
    int numOfRecordFound = 0 ;
    float Key = 8;
    link_list_head *llhead = find(root, Key, true, NULL,&numOfAccess);
    link_list *ll = llhead->head;
    printf("Number of index block access : %d \n",numOfAccess+1);
    if(ll != NULL){
        printf("tconst of record found : \n");
        while(ll!=NULL){
            //hide result too long
            //printf("%s|" , ll->ptr->title);
            numOfRecordFound++;
            ll=ll->nextItem;
        }
    }
    printf("Number of Record Found : %d \n",numOfRecordFound);
    int blockaccess = calculateblockAccess(llhead,head,false);
    printf("Number of data block access : %d \n",blockaccess);
    printf("---------------------------------\n");


    printf("\nRange of Key Search Statistics\n");
    printf("---------------------------------\n");
    //Retrieve for range of key
    //Hard Code searching for range value -> can change to let user type the value
    numOfAccess = 0;
    blockaccess = 0;
    int numOfHead = 0;
    float returned_keys[50000];
	void *returned_pointers[50000];
	float sKey = 7;
	float eKey = 9;
	numOfRecordFound = find_range(root, sKey, eKey, false,returned_keys, returned_pointers,&numOfAccess,&numOfHead);
    printf("Number of index block access : %d \n",numOfAccess+1);
    printf("Number of Record Found : %d \n",numOfRecordFound);
    blockaccess = calculateblockAccessRange(sKey, eKey,head,false,numOfHead);
    printf("Number of data block access : %d \n",blockaccess);
    printf("---------------------------------\n");
    fclose(fp);
    //if (line)
        //free(line);
    //exit(EXIT_SUCCESS);

    //find_and_print(root, 4.1, false);
    //

}


