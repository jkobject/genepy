import gspread
from oauth2client.service_account import ServiceAccountCredentials

scope = ["https://spreadsheets.google.com/feeds", 'https://www.googleapis.com/auth/spreadsheets',
         "https://www.googleapis.com/auth/drive.file", "https://www.googleapis.com/auth/drive"]

credentials = ServiceAccountCredentials.from_json_keyfile_name('~/.client_secret.json', scope)
client = gspread.authorize(credentials)

spreadsheet = client.open('https://docs.google.com/spreadsheets/d/1XkZypRuOEXzNLxVk9EOHeWRE98Z8_DBvL4PovyM01FE')

with open(file, 'r') as file_obj:
  content = file_obj.read()
  client.import_csv(spreadsheet.id, data=content)
