
#include "DataStructures.h"
#include "Constants.h"

lst_item* createListItem(int value){
    struct lst_item* new_node = (struct lst_item*) malloc(sizeof(struct lst_item));
    new_node->next  = NULL;
    new_node->value = value;
    return new_node;
}

lst_item** iniAList(int size){
    lst_item** aList = ( struct lst_item**)malloc(size * sizeof(struct lst_item*));
    for(int i=0; i<size; i++){
        //cria a cabeça de cada lista
        aList[i] = NULL;
    }
    return aList;
}

void free_list(lst_item *head) {
    lst_item *prev = head;
    lst_item *cur = head;
    while(cur) {
        prev = cur;
        cur = prev->next;
        free(prev);
    }       
}

void freeAList(lst_item** aList, int size){
    for(int i=0; i<size; i++){
        free_list(aList[i]);
    }  
    free(aList);
}

/**adiciona item na lista de maneira ordenada*/
void addItem(lst_item** head, int value){
    lst_item* new_node = createListItem(value);
   struct lst_item* current;
    /* Special case for the head end */
    if (*head == NULL || (*head)->value > new_node->value)
    {
        new_node->next = *head;
        *head = new_node;
    }else {
        /* Locate the node before the point of insertion */
        current = *head;
        while (current->next!=NULL &&
               current->next->value < new_node->value)
        {
            //se ja existe na lista, nao adiciiona
            if(current->value==new_node->value){
                free(new_node);
                return;
            }
            current = current->next;
        }
        //so adiciona o novo item na lista se o valor for difrente
        if(current->value!=new_node->value){
            new_node->next = current->next;
            current->next = new_node;
        }else{
            free(new_node);
        }        
    }
}

void printList(lst_item* head){
    printf(":: ");
    struct lst_item *temp = head;
    while(temp != NULL)
    {
        printf("%d  ", temp->value);
        temp = temp->next;
    }
    printf(" ::\n");
}
