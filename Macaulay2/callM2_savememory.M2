callM2 = method();
callM2(String,Boolean) := String=> (command,bQuiet) -> (
     M2io := openInOut concatenate("!M2 " | " --silent "); 
     M2io  << command;
     M2io  << closeOut;
);

innerLoop = (a,b)->(  a+b);

saveResult=(result)->
(
 g = openOut "savedResult.M2";
 g << "result =" | toExternalString result << ";";
 g << close;
);

end


a:=5; b:=7;
command := "load(\"callM2.M2\"); result := innerLoop(" | toExternalString a |","|  toExternalString b |"); saveResult(result) ;" 
callM2( command,false )
load "savedResult.M2"
result
