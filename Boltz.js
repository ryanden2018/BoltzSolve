function hermiteEval(x,n)
{
    var a = Math.exp(-0.5*x*x) * Math.pow(3.141592653589793238462,-0.25);
    var b = Math.sqrt(2.0) * x * Math.exp(-0.5*x*x) * Math.pow(3.141592653589793238462,-0.25);
    for(var i = 0; i < n; i++)
    {
        var c =  Math.sqrt(2.0) * Math.pow((i+2.0),-0.5) * x * b - Math.pow(i+1.0,0.5) * Math.pow(i+2.0,-0.5) * a;
        a = b;
        b = c;
    }
    return a;
}

function hermiteDerivEval(x,n)
{
    var a = Math.exp(-0.5*x*x) * Math.pow(3.141592653589793238462,-0.25);
    var b = Math.sqrt(2.0) * x * Math.exp(-0.5*x*x) * Math.pow(3.141592653589793238462,-0.25);
    var c = x * b - Math.pow(2.0,-0.5) * a;
    for(var i = 0; i < n-1; i++)
    {
        a = b;
        b = c;
        c = Math.pow(2.0,0.5) * Math.pow(i+3.0,-0.5) * x * b - Math.pow(i+2.0,0.5) * Math.pow(i+3.0,-0.5) * a;
    }
    return (n <= 0
    ?
    (-1.0) * x * Math.exp(-0.5 * x * x) * Math.pow(3.141592653589793238462,-0.25)
    :
    Math.pow(n,0.5)*Math.pow(2.0,-0.5)*a-Math.pow(n+1.0,0.5)*Math.pow(2.0,-0.5)*c);
}