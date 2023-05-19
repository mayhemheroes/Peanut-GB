gb-fuzz: gb-fuzz.c
	clang -g -fsanitize=fuzzer,undefined,address -fno-sanitize-recover=undefined -o gb-fuzz gb-fuzz.c