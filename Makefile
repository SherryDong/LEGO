main:
	gcc -o bin/LEGO src/LEGO.c -lm 	
	gcc -o bin/LEGO_pre src/LEGO_pre.c -lm 	
	gcc -o bin/LEGO_mul src/LEGO_mul.c -lm 	
	g++ -o bin/extract_info src/extract_info.cpp -lm	
	g++ -o bin/extract_info_gs src/extract_info_gs.cpp -lm	
