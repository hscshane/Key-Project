from django.shortcuts import render
from django.http import HttpResponse
from .models import SWCX
import requests
from subprocess import run,PIPE
import sys

external_dir = 'D:\\OneDrive - University of Miami\\work\\COSMOS\\swcx\\Sicong\\swcx\\'

def home(request):    
    return render(request, "home.html")

# def result(request):
#     data=requests.get("https://reqres.in//api//users")
#     data = data.text
#     return render(request, 'result.html', {
#         'data': data,
#     })

def external(request):
    coords = request.POST.get('Coordinates')
    coordstype = request.POST.get('Ctype')
    date = request.POST.get('Date')
    #CSRatio = request.POST.get('CSRatio')
    out = run([sys.executable, external_dir+'swcx.py', coords, coordstype, date], shell=False, stdout=PIPE, text=True)

    return render(request, 'home.html', {
        'data1': out.stdout
    })