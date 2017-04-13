CC=gcc

main:
	$(CC) zeta_serial_main.c -o zeta_serial

utest:
	$(CC) unit_test_zeta.c -o ut_zeta
	./ut_zeta

vtest:
	$(CC) ver_test_zeta.c -o vtest_zeta
	./vtest_zeta.exe
