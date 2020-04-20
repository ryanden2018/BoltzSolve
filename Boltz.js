var gaussQuad =[
    [ 0.1392518728556320, -0.0697392733197222],
    [ 0.1392518728556320, 0.0697392733197222],
    [ 0.1365414983460152, -0.2078604266882213],
    [ 0.1365414983460152, 0.2078604266882213],
    [ 0.1311735047870624, -0.3419358208920842],
    [ 0.1311735047870624, 0.3419358208920842],
    [ 0.1232523768105124, -0.4693558379867570],
    [ 0.1232523768105124,  0.4693558379867570],
    [ 0.1129322960805392, -0.5876404035069116],
    [ 0.1129322960805392, 0.5876404035069116],
    [ 0.1004141444428810, -0.6944872631866827],
    [ 0.1004141444428810, 0.6944872631866827],
    [ 0.0859416062170677, -0.7878168059792081],
    [ 0.0859416062170677, 0.7878168059792081],
    [ 0.0697964684245205, -0.8658125777203002],
    [ 0.0697964684245205, 0.8658125777203002],
    [ 0.0522933351526833, -0.9269567721871740],
    [ 0.0522933351526833, 0.9269567721871740],
    [ 0.0337749015848142, -0.9700604978354287],
    [ 0.0337749015848142, 0.9700604978354287],
    [ 0.0146279952982722, -0.9942945854823992],
    [ 0.0146279952982722, 0.9942945854823992]
];

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

function boltzQ()
{
    var ix = this.thread.x%10;
    var iy = (this.thread.x-ix)/10;
    var jx = this.thread.y%10;
    var jy = (this.thread.y-jx)/10;
    var kx = this.thread.z%10;
    var ky = (this.thread.z-kx)/10;
    return 1.85;
}