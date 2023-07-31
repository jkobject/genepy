import gspread
from oauth2client.service_account import ServiceAccountCredentials

scope = ["https://spreadsheets.google.com/feeds", 'https://www.googleapis.com/auth/spreadsheets',
         "https://www.googleapis.com/auth/drive.file", "https://www.googleapis.com/auth/drive"]


def dfToSheet(df, sheetid, secret='~/.credentials.json'):
  credentials = ServiceAccountCredentials.from_json_keyfile_name(secret, scope)
  client = gspread.authorize(credentials)
  spreadsheet = client.open(sheetid)
  df.to_csv('/tmp/sheet.csv')
  with open("/tmp/sheet.csv", 'r') as file_obj:
    content = file_obj.read()
    client.import_csv(spreadsheet.id, data=content)
