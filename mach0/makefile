CC=gcc

main:
	$(CC) machin_serial_main.c -o machin_serial

utest:
	$(CC) unit_test_machin.c -o ut_machin
	./ut_machin
vtest:
	$(CC) ver_test_machin.c -o vtest_machin
	./vtest_machin