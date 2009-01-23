CREATE FUNCTION plpgsql_call_handler() RETURNS OPAQUE AS '/usr/lib64/pgsql/plpgsql.so'  LANGUAGE 'C';
 CREATE LANGUAGE 'plpgsql' HANDLER plpgsql_call_handler
                            LANCOMPILER 'PL/pgSQL';
