build:
	@gcc -c qsufsort.c -o qsufsort.o
	@gcc -c suftest.c -o suftest.o
	@gcc -o suftest suftest.o qsufsort.o

clean:
	@find . -name "*.o" -delete
	@find . -name "suftest" -delete

quality:
	gcc -Wall suftest.c
	gcc -Wall qsufsort.c

test_bible:
	@./suftest corpus/bible.txt

test_world:
	@./suftest corpus/world192.txt
